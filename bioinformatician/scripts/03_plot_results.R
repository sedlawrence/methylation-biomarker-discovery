#!/usr/bin/env Rscript

# GSE59685 plotting and annotation script
#
# Inputs:
# - results/differential_analysis/limma_dmps.tsv
# - results/differential_analysis/samples_used.tsv
# - data/GSE59685_betas.csv.gz
# - data/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz
#
# Outputs:
# - annotated_results.tsv
# - significant_hits.tsv
# - top_dmp_table.tsv
# - volcano_plot.png
# - qq_plot.png
# - pvalue_histogram.png
# - top_cpg_beta_distribution.png
# - cell_proportions_by_case.png
# - heatmap_top_50_dmps.png
# - manhattan_plot.png
# - regional_plot_top_locus.png
# - regional_plot_known_ad_locus.png (only if a known AD locus is present)
# - annotation_summary.png
#
# The manifest already contains the main probe-level annotations we need:
# gene names, RefGene region, CpG island relation, enhancer flag, and regulatory
# feature columns. That lets us keep the script self-contained.

suppressWarnings(suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
}))

library(optparse)

option_list <- list(
  make_option(c("--results"), type = "character", default = "results/differential_analysis/limma_dmps.tsv",
              help = "Differential methylation results table [default %default]"),
  make_option(c("--samples"), type = "character", default = "results/differential_analysis/samples_used.tsv",
              help = "Sample metadata used in the model [default %default]"),
  make_option(c("--betas"), type = "character", default = "data/GSE59685_betas.csv.gz",
              help = "Processed beta matrix from GEO [default %default]"),
  make_option(c("--manifest"), type = "character", default = "data/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz",
              help = "Illumina 450K manifest [default %default]"),
  make_option(c("--out_dir"), type = "character", default = "results/differential_analysis/plots",
              help = "Output directory [default %default]"),
  make_option(c("--fdr_threshold"), type = "double", default = 0.05,
              help = "Strict hit threshold on adjusted p-value [default %default]"),
  make_option(c("--delta_beta_threshold"), type = "double", default = 0.05,
              help = "Strict hit threshold on absolute beta difference [default %default]"),
  make_option(c("--top_n"), type = "integer", default = 25L,
              help = "How many ranked probes to keep in the top DMP table [default %default]"),
  make_option(c("--label_n"), type = "integer", default = 10L,
              help = "How many probes to label on the volcano plot [default %default]"),
  make_option(c("--top_cpg_n"), type = "integer", default = 6L,
              help = "How many probes to show in the beta distribution plot [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(user_lib)) {
  user_lib <- file.path(
    path.expand("~"),
    "Library/R",
    paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = "."),
    "library"
  )
}
if (file.access(user_lib, 2) != 0) {
  user_lib <- file.path(getwd(), ".Rlibrary")
}
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_lib, .libPaths())))

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  install.packages(missing, repos = "https://cloud.r-project.org", lib = user_lib, dependencies = TRUE)
}

install_if_missing(c("data.table"))

suppressPackageStartupMessages({
  library(data.table)
})

read_gz_header <- function(path, n = 10L) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

find_header_line <- function(path, prefix = "IlmnID,") {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  line_no <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    line_no <- line_no + 1L
    if (startsWith(line, prefix)) return(line_no)
  }
  stop("Could not find header line starting with ", prefix, " in: ", path)
}

trim_or_na <- function(x) {
  x <- trimws(x)
  x[!nzchar(x)] <- NA_character_
  x
}

first_token <- function(x) {
  x <- trim_or_na(x)
  if (is.na(x)) return(NA_character_)
  tokens <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) NA_character_ else tokens[1L]
}

classify_gene_region <- function(x) {
  x <- trim_or_na(x)
  if (is.na(x)) return("Intergenic")
  tokens <- unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))
  tokens <- tokens[nzchar(tokens)]
  if (any(tokens %in% c("TSS1500", "TSS200"))) return("Promoter")
  if (any(tokens == "5'UTR")) return("5' UTR")
  if (any(tokens == "1stExon")) return("First exon")
  if (any(tokens == "Body")) return("Gene body")
  if (any(tokens == "3'UTR")) return("3' UTR")
  "Intergenic"
}

classify_cpg_context <- function(x) {
  x <- trim_or_na(x)
  if (is.na(x)) return("Unannotated")
  if (x == "Island") return("CpG island")
  if (x %in% c("N_Shore", "S_Shore")) return("CpG shore")
  if (x %in% c("N_Shelf", "S_Shelf")) return("CpG shelf")
  if (x == "OpenSea") return("Open sea")
  x
}

classify_regulatory <- function(enhancer_flag, regulatory_group) {
  regulatory_group <- trim_or_na(regulatory_group)
  if (isTRUE(enhancer_flag)) return("Enhancer")
  if (!is.na(regulatory_group)) return(regulatory_group)
  "Unannotated"
}

collapse_nonempty <- function(...) {
  vals <- c(...)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0L) return(NA_character_)
  paste(unique(vals), collapse = " | ")
}

short_label <- function(x, n = 80L) {
  x <- as.character(x)
  x <- trimws(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 3L), "..."), x)
}

matches_gene_set <- function(gene_field, gene_set) {
  if (length(gene_field) == 0L || is.na(gene_field) || !nzchar(gene_field)) return(FALSE)
  tokens <- trimws(toupper(unlist(strsplit(as.character(gene_field), ";", fixed = TRUE))))
  tokens <- tokens[nzchar(tokens)]
  any(tokens %in% toupper(gene_set))
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", x, perl = TRUE)
}

read_manifest_annotations <- function(path) {
  header_line <- find_header_line(path)
  cols <- c(
    "IlmnID", "CHR", "MAPINFO",
    "UCSC_RefGene_Name", "UCSC_RefGene_Group",
    "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island",
    "Enhancer", "Regulatory_Feature_Name", "Regulatory_Feature_Group"
  )
  man <- suppressWarnings(fread(
    path,
    skip = header_line - 1L,
    header = TRUE,
    select = cols,
    showProgress = FALSE
  ))
  setnames(man, "IlmnID", "probe_id")
  man[, gene_name := vapply(UCSC_RefGene_Name, first_token, character(1))]
  man[, gene_region := vapply(UCSC_RefGene_Group, classify_gene_region, character(1))]
  man[, cpg_context := vapply(Relation_to_UCSC_CpG_Island, classify_cpg_context, character(1))]
  man[, enhancer_flag := as.logical(Enhancer)]
  man[, regulatory_context := mapply(classify_regulatory, enhancer_flag, Regulatory_Feature_Group, USE.NAMES = FALSE)]
  man[, annotation_summary := mapply(
    collapse_nonempty,
    gene_region,
    cpg_context,
    regulatory_context,
    USE.NAMES = FALSE
  )]
  man[, c(
    "gene_name_all",
    "gene_region_raw",
    "cpg_context_raw"
  ) := .(UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_UCSC_CpG_Island)]
  man[]
}

read_betas_for_probes <- function(path, probe_ids, barcodes_keep) {
  hdr <- read_gz_header(path, n = 7L)
  if (length(hdr) < 7L) stop("Unexpected beta file header: ", path)
  sample_ids <- strsplit(sub("^,", "", hdr[5]), ",", fixed = TRUE)[[1]]
  idx <- match(barcodes_keep, sample_ids)
  if (anyNA(idx)) {
    stop("Some requested barcodes are missing from the beta file: ",
         paste(barcodes_keep[is.na(idx)], collapse = ", "))
  }

  probe_ids <- unique(probe_ids)
  probe_ids <- probe_ids[nzchar(probe_ids)]
  if (length(probe_ids) == 0L) stop("No probe IDs requested from beta matrix.")

  gzcat <- Sys.which("gzcat")
  if (!nzchar(gzcat)) gzcat <- Sys.which("zcat")
  grep_bin <- Sys.which("grep")
  tail_bin <- Sys.which("tail")

  if (nzchar(gzcat) && nzchar(grep_bin) && nzchar(tail_bin)) {
    pattern <- paste0("^(", paste(escape_regex(probe_ids), collapse = "|"), "),")
    cmd <- sprintf(
      "%s %s | %s -n +8 | %s -E %s",
      shQuote(gzcat),
      shQuote(path),
      shQuote(tail_bin),
      shQuote(grep_bin),
      shQuote(pattern)
    )
    dt <- fread(cmd = cmd, sep = ",", header = FALSE, showProgress = FALSE)
  } else {
    dt <- fread(path, sep = ",", skip = 7L, header = FALSE, showProgress = FALSE)
    dt <- dt[V1 %in% probe_ids]
  }

  if (nrow(dt) == 0L) stop("Could not find any requested probes in: ", path)
  setnames(dt, c("probe_id", sample_ids))
  dt <- dt[probe_id %in% probe_ids, c("probe_id", barcodes_keep), with = FALSE]
  setnames(dt, c("probe_id", barcodes_keep))
  dt <- dt[match(probe_ids, probe_id)]
  dt <- dt[!is.na(probe_id)]
  dt
}

save_plot <- function(filename, width, height, expr) {
  png(
    filename = file.path(opt$out_dir, filename),
    width = width,
    height = height,
    units = "in",
    res = 300,
    bg = "white"
  )
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

chr_sort_key <- function(chr) {
  chr <- toupper(trimws(as.character(chr)))
  out <- rep(NA_integer_, length(chr))
  num <- suppressWarnings(as.integer(chr))
  out[!is.na(num)] <- num[!is.na(num)]
  out[chr == "X"] <- 23L
  out[chr == "Y"] <- 24L
  out[chr == "XY"] <- 25L
  out[chr %in% c("M", "MT", "CHRM")] <- 26L
  out[is.na(out)] <- 99L
  out
}

make_qq_plot <- function(pvals, filename) {
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  n <- length(pvals)
  if (n == 0L) return(invisible(NULL))
  obs <- -log10(sort(pvals))
  exp <- -log10(ppoints(n))
  lambda <- median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

  save_plot(filename, 7.5, 7, {
    op <- par(mar = c(5, 5, 4, 2) + 0.1)
    on.exit(par(op), add = TRUE)
    plot(exp, obs,
      pch = 16, cex = 0.5, col = "grey40",
      xlab = expression(Expected~(-log[10](P))),
      ylab = expression(Observed~(-log[10](P))),
      main = "QQ plot of raw p-values"
    )
    abline(0, 1, col = "#d73027", lwd = 2)
    mtext(sprintf("lambda = %.3f", lambda), side = 3, adj = 1, line = 0.2, cex = 0.8)
  })
}

make_p_hist <- function(pvals, filename) {
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  save_plot(filename, 7.5, 6.5, {
    op <- par(mar = c(5, 5, 4, 2) + 0.1)
    on.exit(par(op), add = TRUE)
    hist(
      pvals,
      breaks = 50,
      probability = TRUE,
      col = "grey75",
      border = "white",
      main = "Raw p-value histogram",
      xlab = "Raw p-value"
    )
    abline(h = 1, col = "#d73027", lwd = 2, lty = 2)
  })
}

make_manhattan_plot <- function(annotated_dt, top_table_dt, filename) {
  man <- annotated_dt[!is.na(CHR) & !is.na(MAPINFO) & is.finite(P.Value)]
  if (nrow(man) == 0L) return(invisible(NULL))
  man[, chr_chr := toupper(trimws(as.character(CHR)))]
  man[, chr_key := chr_sort_key(chr_chr)]
  setorder(man, chr_key, MAPINFO)

  chr_levels <- unique(man$chr_chr)
  chr_summary <- man[, .(
    chr_key = min(chr_key),
    min_pos = min(MAPINFO),
    max_pos = max(MAPINFO)
  ), by = chr_chr][order(chr_key)]
  chr_summary[, offset := c(0, cumsum(head(max_pos, -1) + 1e6))]
  man <- merge(man, chr_summary[, .(chr_chr, offset)], by = "chr_chr", all.x = TRUE, sort = FALSE)
  man[, x := MAPINFO + offset]
  man[, y := neglog10_p]
  man[, point_col := ifelse(chr_key %% 2 == 0, "grey35", "grey70")]

  top_lab <- top_table_dt[1:min(10L, nrow(top_table_dt))]
  top_lab <- man[probe_id %in% top_lab$probe_id]
  top_lab[, plot_label := ifelse(
    !is.na(gene_name) & nzchar(gene_name),
    paste0(probe_id, "\n", gene_name),
    probe_id
  )]

  save_plot(filename, 10.5, 6.5, {
    op <- par(mar = c(5, 5, 4, 2) + 0.1, xaxs = "i")
    on.exit(par(op), add = TRUE)
    plot(
      man$x, man$y,
      pch = 16, cex = 0.35,
      col = man$point_col,
      xlab = "Genomic position",
      ylab = expression(-log[10](P.value)),
      main = "Manhattan plot"
    )
    if (nrow(top_lab) > 0L) {
      points(top_lab$x, top_lab$y, pch = 16, cex = 0.55, col = "#d73027")
      text(top_lab$x, top_lab$y, labels = top_lab$plot_label, pos = 3, cex = 0.55, xpd = NA)
    }
    chr_mid <- chr_summary[, .(mid = (min_pos + max_pos) / 2 + offset, chr_chr)]
    axis(1, at = chr_mid$mid, labels = chr_mid$chr_chr, las = 2, cex.axis = 0.75)
    legend("topright", legend = c("Background", "Top loci"), col = c("grey35", "#d73027"),
           pch = 16, bty = "n", cex = 0.85)
  })
}

make_regional_plot <- function(annotated_dt, center_probe, filename, window_bp = 250000L, label_n = 8L) {
  center <- annotated_dt[probe_id == center_probe]
  if (nrow(center) == 0L || is.na(center$CHR[1]) || is.na(center$MAPINFO[1])) return(invisible(FALSE))
  chr_chr <- toupper(trimws(as.character(center$CHR[1])))
  center_pos <- as.numeric(center$MAPINFO[1])
  region <- annotated_dt[
    toupper(trimws(as.character(CHR))) == chr_chr &
      MAPINFO >= (center_pos - window_bp) &
      MAPINFO <= (center_pos + window_bp)
  ]
  if (nrow(region) == 0L) return(invisible(FALSE))
  setorder(region, MAPINFO)
  region[, y := neglog10_p]
  region[, is_center := probe_id == center_probe]
  region[, plot_col := fifelse(is_center, "#d73027", "grey50")]

  label_dt <- region[order(P.Value)][1:min(label_n, nrow(region))]
  label_dt[, plot_label := ifelse(
    !is.na(gene_name) & nzchar(gene_name),
    paste0(probe_id, "\n", gene_name),
    probe_id
  )]

  title_txt <- if (!is.na(center$gene_name[1]) && nzchar(center$gene_name[1])) {
    paste0("Regional plot: ", center_probe, " | ", center$gene_name[1])
  } else {
    paste0("Regional plot: ", center_probe)
  }

  save_plot(filename, 10.5, 6, {
    op <- par(mar = c(5, 5, 4, 2) + 0.1, xaxs = "i")
    on.exit(par(op), add = TRUE)
    plot(
      region$MAPINFO, region$y,
      pch = 16, cex = 0.6,
      col = region$plot_col,
      xlab = paste0("Position on chr ", chr_chr, " (bp)"),
      ylab = expression(-log[10](P.value)),
      main = title_txt
    )
    abline(v = center_pos, col = "#d73027", lty = 2)
    if (nrow(label_dt) > 0L) {
      text(label_dt$MAPINFO, label_dt$y, labels = label_dt$plot_label, pos = 3, cex = 0.6, xpd = NA)
    }
    legend("topright", legend = c("Center probe", "Nearby probes"), col = c("#d73027", "grey50"),
           pch = 16, bty = "n", cex = 0.85)
  })

  invisible(TRUE)
}

message("[1/6] Loading results and sample metadata...")
res <- fread(opt$results)
samples <- fread(opt$samples)

required_res <- c("probe_id", "P.Value", "adj.P.Val", "delta_beta_model", "delta_beta_mean", "se_M")
missing_res <- setdiff(required_res, names(res))
if (length(missing_res) > 0L) {
  stop("Results table is missing required columns: ", paste(missing_res, collapse = ", "))
}
if (!"case" %in% names(samples)) {
  stop("Sample table is missing 'case'.")
}

message("[2/6] Loading manifest annotations and merging to results...")
manifest <- read_manifest_annotations(opt$manifest)
setkey(manifest, probe_id)
setkey(res, probe_id)

annotated <- manifest[res]
annotated[, hit := adj.P.Val < opt$fdr_threshold & abs(delta_beta_model) >= opt$delta_beta_threshold]
annotated[, neglog10_p := -log10(pmax(P.Value, .Machine$double.xmin))]
annotated[, neglog10_adj_p := -log10(pmax(adj.P.Val, .Machine$double.xmin))]
annotated <- annotated[order(adj.P.Val, -abs(delta_beta_model), P.Value)]

strict_hits <- annotated[hit == TRUE]
if (nrow(strict_hits) == 0L) {
  top_table <- annotated[1:min(opt$top_n, nrow(annotated))]
  top_table[, selected_reason := "top_ranked_fallback"]
  message("  No strict hits passed adj.P.Val < ", opt$fdr_threshold,
          " and |delta_beta_model| >= ", opt$delta_beta_threshold,
          ". Using top-ranked probes for the summary table and plots.")
} else {
  top_table <- strict_hits[1:min(opt$top_n, nrow(strict_hits))]
  top_table[, selected_reason := "strict_hit"]
  message("  Strict hits: ", nrow(strict_hits))
}

annotated[, selected_reason := fifelse(hit, "strict_hit", "not_selected")]

message("[3/6] Writing annotation tables...")
fwrite(annotated, file.path(opt$out_dir, "annotated_results.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(strict_hits, file.path(opt$out_dir, "significant_hits.tsv"), sep = "\t", quote = FALSE, na = "NA")

top_table_out <- top_table[, .(
  probe_id,
  gene_name,
  gene_name_all,
  gene_region,
  gene_region_raw,
  cpg_context,
  cpg_context_raw,
  enhancer_flag,
  regulatory_context,
  Regulatory_Feature_Name,
  Regulatory_Feature_Group,
  CHR,
  MAPINFO,
  logFC_M,
  se_M,
  t,
  P.Value,
  adj.P.Val,
  delta_beta_mean,
  delta_beta_model,
  se_delta_beta_model,
  mean_beta_ctrl,
  mean_beta_case,
  hit,
  selected_reason,
  annotation_summary
)]
fwrite(top_table_out, file.path(opt$out_dir, "top_dmp_table.tsv"), sep = "\t", quote = FALSE, na = "NA")

analysis_samples <- samples[!is.na(case)]
if (nrow(analysis_samples) == 0L) {
  stop("No case/control samples available for plotting.")
}

message("[4/8] Making raw p-value diagnostics...")
make_qq_plot(res$P.Value, "qq_plot.png")
make_p_hist(res$P.Value, "pvalue_histogram.png")

message("[5/8] Making volcano plot...")
label_rows <- if (nrow(strict_hits) > 0L) {
  strict_hits[1:min(opt$label_n, nrow(strict_hits))]
} else {
  top_table[1:min(opt$label_n, nrow(top_table))]
}
label_rows[, plot_label := ifelse(
  !is.na(gene_name) & nzchar(gene_name),
  paste0(probe_id, "\n", gene_name),
  probe_id
)]

plot_data <- annotated
if (nrow(plot_data) > 100000L) {
  hit_plot <- plot_data[hit == TRUE]
  bg_plot <- plot_data[hit == FALSE]
  if (nrow(bg_plot) > 50000L) {
    set.seed(1)
    bg_plot <- bg_plot[sample.int(nrow(bg_plot), 50000L)]
  }
  plot_data <- rbind(hit_plot, bg_plot, fill = TRUE)
}
plot_data <- unique(rbind(plot_data, label_rows, fill = TRUE), by = "probe_id")

save_plot("volcano_plot.png", 8.5, 6.5, {
  op <- par(mar = c(5, 5, 4, 2) + 0.1, xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot(
    plot_data$delta_beta_model,
    plot_data$neglog10_adj_p,
    pch = 16,
    cex = 0.45,
    col = ifelse(plot_data$hit, "#d73027", "grey70"),
    xlab = expression(Delta * beta[model]),
    ylab = expression(-log[10](adj.P.Val)),
    main = "GSE59685 blood DMP volcano plot"
  )
  abline(v = c(-opt$delta_beta_threshold, opt$delta_beta_threshold), lty = 2, col = "grey45")
  abline(h = -log10(opt$fdr_threshold), lty = 2, col = "grey45")
  legend("topright", legend = c("Not hit", "Strict hit"), col = c("grey70", "#d73027"),
         pch = 16, bty = "n", cex = 0.9)
  if (nrow(label_rows) > 0L) {
    text(
      label_rows$delta_beta_model,
      label_rows$neglog10_adj_p,
      labels = label_rows$plot_label,
      pos = 4,
      cex = 0.6,
      xpd = NA
    )
  }
  if (nrow(strict_hits) == 0L) {
    mtext("No strict hits; labels show top-ranked probes", side = 3, line = 0.2, cex = 0.8)
  }
})

message("[6/8] Making beta boxplots, cell composition plots, and heatmap...")

beta_plot_probes <- top_table$probe_id[seq_len(min(opt$top_cpg_n, nrow(top_table)))]
beta_dt <- read_betas_for_probes(opt$betas, beta_plot_probes, analysis_samples$barcode)
beta_long <- melt(beta_dt, id.vars = "probe_id", variable.name = "barcode", value.name = "beta")
beta_long <- merge(beta_long, analysis_samples[, .(barcode, case)], by = "barcode", all.x = TRUE, sort = FALSE)
beta_long[, case_label := factor(case, levels = c(0, 1), labels = c("Control", "Case"))]
beta_long <- merge(
  beta_long,
  top_table[, .(probe_id, gene_name, gene_region, cpg_context, annotation_summary)],
  by = "probe_id",
  all.x = TRUE,
  sort = FALSE
)

save_plot("top_cpg_beta_distribution.png", 11, 8.5, {
  nplots <- length(beta_plot_probes)
  ncol <- 2L
  nrow <- ceiling(nplots / ncol)
  op <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 2, 0))
  on.exit(par(op), add = TRUE)
  for (pid in beta_plot_probes) {
    sub <- beta_long[probe_id == pid]
    sub[, case_label := factor(case_label, levels = c("Control", "Case"))]
    x <- split(sub$beta, sub$case_label)
    main_label <- short_label(
      if (all(is.na(sub$gene_name))) {
        paste(pid, sub$gene_region[1], sub$cpg_context[1], sep = " | ")
      } else {
        paste(pid, sub$gene_name[1], sub$gene_region[1], sub$cpg_context[1], sep = " | ")
      },
      n = 70L
    )
    boxplot(
      x,
      ylim = c(0, 1),
      col = c("#4575b4", "#d73027"),
      border = "grey30",
      ylab = "Beta",
      xlab = "",
      main = main_label
    )
    stripchart(
      x,
      vertical = TRUE,
      method = "jitter",
      pch = 16,
      col = c(adjustcolor("#4575b4", alpha.f = 0.45), adjustcolor("#d73027", alpha.f = 0.45)),
      add = TRUE
    )
    abline(h = 0.5, lty = 3, col = "grey80")
  }
  mtext("Beta distributions for top ranked probes", outer = TRUE, cex = 1.2, font = 2)
  mtext("Each panel shows raw beta values in case/control groups", outer = TRUE, line = -1, cex = 0.8)
})

cell_cols <- intersect(c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"), names(analysis_samples))
if (length(cell_cols) > 0L) {
  cell_long <- melt(
    analysis_samples,
    id.vars = c("barcode", "case"),
    measure.vars = cell_cols,
    variable.name = "cell_type",
    value.name = "proportion"
  )
  cell_long[, case_label := factor(case, levels = c(0, 1), labels = c("Control", "Case"))]
  save_plot("cell_proportions_by_case.png", 11, 7.5, {
    nplots <- length(cell_cols)
    ncol <- 3L
    nrow <- ceiling(nplots / ncol)
    op <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 2, 0))
    on.exit(par(op), add = TRUE)
    for (ct in cell_cols) {
      sub <- cell_long[cell_type == ct]
      x <- split(sub$proportion, sub$case_label)
      boxplot(
        x,
        ylim = c(0, max(1, max(sub$proportion, na.rm = TRUE) * 1.1)),
        col = c("#4575b4", "#d73027"),
        border = "grey30",
        ylab = "Estimated proportion",
        xlab = "",
        main = ct
      )
      stripchart(
        x,
        vertical = TRUE,
        method = "jitter",
        pch = 16,
        col = c(adjustcolor("#4575b4", alpha.f = 0.45), adjustcolor("#d73027", alpha.f = 0.45)),
        add = TRUE
      )
    }
    mtext("Estimated blood cell proportions by case/control", outer = TRUE, cex = 1.2, font = 2)
    mtext("Derived from the Houseman-style deconvolution used in the analysis", outer = TRUE, line = -1, cex = 0.8)
  })
}

top_heatmap <- annotated[1:min(50L, nrow(annotated))]
top_heatmap <- top_heatmap[!is.na(probe_id)]
heatmap_dt <- read_betas_for_probes(opt$betas, top_heatmap$probe_id, analysis_samples$barcode)
heatmap_mat <- as.matrix(heatmap_dt[, -1, with = FALSE])
rownames(heatmap_mat) <- heatmap_dt$probe_id
heatmap_order <- match(top_heatmap$probe_id, rownames(heatmap_mat))
heatmap_mat <- heatmap_mat[heatmap_order, , drop = FALSE]
top_heatmap <- top_heatmap[match(rownames(heatmap_mat), probe_id)]
sample_order <- order(factor(analysis_samples$case, levels = c(0, 1)), analysis_samples$barcode)
analysis_samples_heat <- analysis_samples[sample_order]
heatmap_mat <- heatmap_mat[, analysis_samples_heat$barcode, drop = FALSE]
col_side <- c("#4575b4", "#d73027")[as.integer(factor(analysis_samples_heat$case, levels = c(0, 1)))]
heatmap_labels <- ifelse(
  !is.na(top_heatmap$gene_name) & nzchar(top_heatmap$gene_name),
  paste0(top_heatmap$probe_id, " | ", top_heatmap$gene_name),
  top_heatmap$probe_id
)

save_plot("heatmap_top_50_dmps.png", 12, 10, {
  op <- par(mar = c(8, 8, 4, 2) + 0.1)
  on.exit(par(op), add = TRUE)
  heatmap(
    heatmap_mat,
    scale = "row",
    Rowv = TRUE,
    Colv = NA,
    col = colorRampPalette(c("#313695", "white", "#a50026"))(101),
    labRow = short_label(heatmap_labels, 60L),
    labCol = FALSE,
    ColSideColors = col_side,
    margins = c(10, 10),
    main = "Top 50 DMPs (row-scaled beta values)"
  )
  legend(
    "topright",
    legend = c("Control", "Case"),
    fill = c("#4575b4", "#d73027"),
    bty = "n",
    cex = 0.85
  )
})

message("[7/8] Making Manhattan and regional plots...")
make_manhattan_plot(annotated, top_table, "manhattan_plot.png")

if (nrow(top_table) > 0L) {
  make_regional_plot(annotated, top_table$probe_id[1], "regional_plot_top_locus.png")
}

known_ad_genes <- c(
  "ANK1", "ABCA7", "APOE", "BIN1", "CLU", "CR1", "PICALM", "CD33",
  "MS4A6A", "MS4A4A", "HLA-DRA", "HLA-DRB1", "INPP5D", "PLCG2", "SORL1",
  "TREM2", "EPHA1", "FERMT2", "CASS4", "RIN3", "MEF2C", "SLC24A4",
  "SPI1", "ABI3", "PLD3", "PTK2B", "PILRA", "CELF1", "NME8", "FYN"
)
known_ad_hits <- annotated[vapply(gene_name_all, matches_gene_set, logical(1), gene_set = known_ad_genes)]
if (nrow(known_ad_hits) > 0L) {
  known_ad_hits <- known_ad_hits[probe_id != top_table$probe_id[1]]
}
if (nrow(known_ad_hits) > 0L) {
  make_regional_plot(annotated, known_ad_hits$probe_id[1], "regional_plot_known_ad_locus.png")
  message("  Known AD locus regional plot created for: ", known_ad_hits$probe_id[1], " (", known_ad_hits$gene_name[1], ")")
} else {
  message("  No additional known AD locus found among the results; skipping the known-AD regional plot.")
}

message("[8/8] Making annotation summary plot...")
save_plot("annotation_summary.png", 11, 8, {
  op <- par(mfrow = c(2, 2), mar = c(7, 4, 3, 1) + 0.1, oma = c(0, 0, 2, 0))
  on.exit(par(op), add = TRUE)
  plots <- list(
    "Gene region" = top_table$gene_region,
    "CpG context" = top_table$cpg_context,
    "Regulatory" = top_table$regulatory_context,
    "Enhancer flag" = ifelse(top_table$enhancer_flag, "Enhancer", "No enhancer")
  )
  for (nm in names(plots)) {
    tbl <- table(plots[[nm]])
    tbl <- sort(tbl, decreasing = FALSE)
    barplot(
      tbl,
      las = 2,
      col = "steelblue",
      border = "grey30",
      main = nm,
      ylab = "Count"
    )
  }
  mtext("Annotation summary for selected probes", outer = TRUE, cex = 1.2, font = 2)
  mtext("Using manifest metadata already present for GSE59685/GPL13534", outer = TRUE, line = -1, cex = 0.8)
})

message("Wrote tables and plots to: ", opt$out_dir)

sink(file.path(opt$out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
