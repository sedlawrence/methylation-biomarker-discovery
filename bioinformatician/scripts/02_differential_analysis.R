#!/usr/bin/env Rscript

# GSE59685 (whole blood only): QC -> cell composition adjustment -> limma DMPs
#
# Notes on cell composition:
# - The canonical minfi workflow (minfi::estimateCellCounts) requires an RGChannelSet
#   built from IDATs. GSE59685 does NOT provide IDATs; it provides processed betas and
#   GEO "signal intensities" tables.
# - To still use the standard Houseman-style reference (FlowSorted.Blood.450k), we:
#     1) preprocess the reference (FlowSorted.Blood.450k) with preprocessQuantile
#     2) select cell-type discriminating probes via minfi:::pickCompProbes
#     3) project *our* beta matrix onto those reference coefficients via
#        minfi:::projectCellType
#   This matches the deconvolution model used inside estimateCellCounts, but cannot
#   reproduce its joint (user+reference) RGChannelSet processing/normalization step.
#
# Outputs (default under results/differential_analysis/):
# - samples_used.tsv
# - cell_proportions.tsv
# - limma_dmps.tsv
# - sessionInfo.txt

suppressWarnings(suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse", repos = "https://cloud.r-project.org")
}))

library(optparse)

option_list <- list(
  make_option(c("--sample_qc"), type = "character", default = "results/gse59685_qc/sample_qc.tsv.gz",
              help = "Sample QC table produced by step 1 [default %default]"),
  make_option(c("--probe_qc"), type = "character", default = "results/gse59685_qc/probe_qc.tsv.gz",
              help = "Probe QC table produced by step 1 [default %default]"),
  make_option(c("--betas"), type = "character", default = "data/GSE59685_betas.csv.gz",
              help = "Processed betas matrix from GEO [default %default]"),
  make_option(c("--out_dir"), type = "character", default = "results/differential_analysis",
              help = "Output directory [default %default]"),
  make_option(c("--braak_control_max"), type = "integer", default = 2L,
              help = "Controls are Braak <= this value [default %default]"),
  make_option(c("--braak_case_min"), type = "integer", default = 3L,
              help = "Cases are Braak >= this value [default %default]"),
  make_option(c("--drop_celltype"), type = "character", default = "Gran",
              help = "Cell type covariate to drop from design to avoid collinearity [default %default]"),
  make_option(c("--no_install"), action = "store_true", default = FALSE,
              help = "Do not install missing packages (assume environment is prepared) [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))

  if (isTRUE(opt$no_install)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nRe-run without --no_install to auto-install.")
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  options(repos = BiocManager::repositories())

  if (bioc) {
    message("Bioconductor repositories in use:")
    print(getOption("repos"))
    BiocManager::install(missing, ask = FALSE, update = TRUE, dependencies = TRUE)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

# CRAN
install_if_missing(c("data.table"))
install_if_missing(c("quadprog"))

# Bioconductor
install_if_missing(c("limma", "minfi", "FlowSorted.Blood.450k", "SummarizedExperiment"), bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(minfi)
  library(FlowSorted.Blood.450k)
})

as_lgl <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  x %in% c("true", "t", "1", "yes", "y")
}

read_gse59685_betas_subset <- function(path, barcodes_keep) {
  # File structure:
  #   lines 1-3: comments
  #   line 4: blank
  #   line 5: ",<barcode1>,<barcode2>,..."
  #   line 6: ",<new GSM1>,..."
  #   line 7: ",<old GSM1>,..."
  #   line 8+: "<cg_id>,<beta1>,<beta2>,..."
  con <- gzfile(path, "rt")
  hdr <- readLines(con, n = 7, warn = FALSE)
  close(con)
  if (length(hdr) < 7) stop("Unexpected betas file header; expected at least 7 lines: ", path)

  barcode_line <- hdr[5]
  barcode_line <- sub("^\uFEFF", "", barcode_line) # strip UTF-8 BOM if present
  barcode_line <- sub("^,", "", barcode_line)
  barcodes_all <- strsplit(barcode_line, ",", fixed = TRUE)[[1]]
  if (length(barcodes_all) < 10) stop("Failed to parse sample barcodes from betas file header.")

  idx <- match(barcodes_keep, barcodes_all)
  if (any(is.na(idx))) {
    missing <- barcodes_keep[is.na(idx)]
    stop("Some requested barcodes are missing from betas file: ", paste(missing, collapse = ", "))
  }
  # Column 1 is probe_id, then sample columns start at 2.
  select_cols <- c(1L, idx + 1L)

  dt <- fread(path, sep = ",", skip = 7, header = FALSE, select = select_cols, showProgress = TRUE)
  setnames(dt, c("probe_id", barcodes_keep))
  dt
}

beta_from_M <- function(M) {
  # beta = 2^M / (1 + 2^M)
  t <- 2^M
  t / (1 + t)
}

estimate_cell_props_houseman <- function(beta_mat,
                                         cell_types = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                                         probe_select = "both") {
  data("FlowSorted.Blood.450k", package = "FlowSorted.Blood.450k")
  ref_rg <- FlowSorted.Blood.450k

  # Match estimateCellCounts default for Blood: preprocessQuantile + probeSelect="both"
  ref_mset <- preprocessQuantile(ref_rg)

  compData <- minfi:::pickCompProbes(
    mSet = ref_mset,
    cellTypes = cell_types,
    compositeCellType = "Blood",
    probeSelect = probe_select
  )
  coefs <- compData$coefEsts

  common <- intersect(rownames(coefs), rownames(beta_mat))
  if (length(common) < 0.8 * nrow(coefs)) {
    warning("Only ", length(common), "/", nrow(coefs),
            " cell-composition training probes are present after QC filtering.")
  }
  y <- beta_mat[common, , drop = FALSE]
  coefs2 <- coefs[common, , drop = FALSE]

  props <- minfi:::projectCellType(Y = y, coefCellType = coefs2)
  props
}

houseman_reference_coefs <- function(cell_types = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                                     probe_select = "both") {
  data("FlowSorted.Blood.450k", package = "FlowSorted.Blood.450k")
  ref_rg <- FlowSorted.Blood.450k
  ref_mset <- preprocessQuantile(ref_rg)

  compData <- minfi:::pickCompProbes(
    mSet = ref_mset,
    cellTypes = cell_types,
    compositeCellType = "Blood",
    probeSelect = probe_select
  )
  compData$coefEsts
}

message("[1/6] Loading QC tables and selecting whole blood samples/probes...")
sample_qc <- fread(opt$sample_qc)
probe_qc <- fread(opt$probe_qc)

sample_qc[, sample_final_keep := as_lgl(sample_final_keep)]
probe_qc[, final_keep := as_lgl(final_keep)]

blood <- sample_qc[source_tissue == "whole blood" & sample_final_keep == TRUE]
if (nrow(blood) == 0) stop("No whole blood samples found with sample_final_keep=TRUE.")

keep_probes <- probe_qc[final_keep == TRUE, probe_id]
if (length(keep_probes) == 0) stop("No probes found with final_keep=TRUE.")

# Braak -> numeric; "Exclude" becomes NA
blood[, braak_num := suppressWarnings(as.integer(braak_stage))]
blood[, case := fifelse(!is.na(braak_num) & braak_num >= opt$braak_case_min, 1L,
                        fifelse(!is.na(braak_num) & braak_num <= opt$braak_control_max, 0L, NA_integer_))]

analysis_samples <- blood[!is.na(case)]
if (nrow(analysis_samples) < 6) {
  stop("Too few samples after Braak binning (n=", nrow(analysis_samples), "). ",
       "Check --braak_control_max/--braak_case_min and braak_stage values.")
}

message("  Whole blood QC-passed samples: ", nrow(blood))
message("  Analysis samples after Braak binning: ", nrow(analysis_samples),
        " (cases=", sum(analysis_samples$case == 1L), ", controls=", sum(analysis_samples$case == 0L), ")")
message("  QC-passed probes for limma: ", length(keep_probes))

message("[2/6] Reading betas for QC whole blood samples (subset of columns only)...")
barcodes_keep <- blood$barcode
beta_dt <- read_gse59685_betas_subset(opt$betas, barcodes_keep = barcodes_keep)

setDT(beta_dt)
setkey(beta_dt, probe_id)

cell_types <- c("CD8T","CD4T","NK","Bcell","Mono","Gran")

message("[3/6] Estimating blood cell proportions (Houseman-style; FlowSorted.Blood.450k reference)...")
ref_coefs <- houseman_reference_coefs(cell_types = cell_types, probe_select = "both")
training_probes <- rownames(ref_coefs)

# Cell type projection only needs the (small) training probe set
beta_comp_dt <- beta_dt[J(training_probes), nomatch = 0]
beta_comp_mat <- as.matrix(beta_comp_dt[, ..barcodes_keep])
rownames(beta_comp_mat) <- beta_comp_dt$probe_id

common <- intersect(training_probes, rownames(beta_comp_mat))
if (length(common) < 0.8 * length(training_probes)) {
  warning("Only ", length(common), "/", length(training_probes),
          " cell-composition training probes were found in the beta matrix.")
}
cell_props <- minfi:::projectCellType(
  Y = beta_comp_mat[common, , drop = FALSE],
  coefCellType = ref_coefs[common, , drop = FALSE]
)

# Persist cell proportions
cell_props_df <- data.frame(barcode = rownames(cell_props), cell_props, check.names = FALSE)
fwrite(cell_props_df, file = file.path(opt$out_dir, "cell_proportions.tsv"), sep = "\t", quote = FALSE, na = "NA")

message("[4/6] Subsetting to final QC probes for differential methylation...")
probe_keep_dt <- data.table(probe_id = keep_probes)
beta_keep_dt <- beta_dt[probe_keep_dt, nomatch = 0]

beta_mat <- as.matrix(beta_keep_dt[, ..barcodes_keep])
rownames(beta_mat) <- beta_keep_dt$probe_id

# Subset to analysis samples (after Braak binning)
analysis_barcodes <- analysis_samples$barcode
beta_mat <- beta_mat[, analysis_barcodes, drop = FALSE]

message("[5/6] Building design matrix (case + age + sex + cell proportions) and running limma on M-values...")
pheno <- as.data.frame(analysis_samples)
rownames(pheno) <- pheno$barcode

# Merge cell props
cell_props_use <- as.data.frame(cell_props)
cell_props_use$barcode <- rownames(cell_props_use)
pheno <- merge(pheno, cell_props_use, by = "barcode", all.x = TRUE, sort = FALSE)
rownames(pheno) <- pheno$barcode

# Ensure covariate rows match methylation matrix columns (and in the same order)
pheno <- pheno[colnames(beta_mat), , drop = FALSE]

# Clean covariates
pheno$age_blood <- as.numeric(pheno$age_blood)
pheno$sex <- factor(pheno$sex)
pheno$case <- as.integer(pheno$case)

# Drop samples with missing covariates
covars <- c("case", "age_blood", "sex", cell_types)
na_any <- !complete.cases(pheno[, covars])
if (any(na_any)) {
  message("  Dropping ", sum(na_any), " samples with missing covariates.")
  keep_rows <- rownames(pheno)[!na_any]
  pheno <- pheno[keep_rows, , drop = FALSE]
  beta_mat <- beta_mat[, keep_rows, drop = FALSE]
}

# M-values (logit base 2)
eps <- 1e-6
beta_clip <- pmin(pmax(beta_mat, eps), 1 - eps)
M <- log2(beta_clip / (1 - beta_clip))

## Design formula
# Avoid collinearity: drop one cell type (default Gran) since proportions sum to ~1.
drop_cell <- opt$drop_celltype
if (!drop_cell %in% cell_types) {
  stop("--drop_celltype must be one of: ", paste(cell_types, collapse = ", "))
}
cell_covars <- setdiff(cell_types, drop_cell)

design_formula <- as.formula(paste("~ case + age_blood + sex +", paste(cell_covars, collapse = " + ")))
design <- model.matrix(design_formula, data = pheno)
if (qr(design)$rank < ncol(design)) {
  stop("Design matrix is rank-deficient (collinearity). Try dropping a different cell type with --drop_celltype.")
}

fit <- lmFit(M, design)
fit <- eBayes(fit)

coef_name <- "case"
if (!coef_name %in% colnames(design)) stop("Expected coefficient not found in design: ", coef_name)
coef_idx <- which(colnames(design) == coef_name)

tt <- topTable(fit, coef = coef_idx, number = Inf, sort.by = "P", adjust.method = "BH")

# Standard error for the coefficient on the (moderated) M-value scale
se_M <- fit$stdev.unscaled[, coef_idx] * sqrt(fit$s2.post)
names(se_M) <- rownames(fit$coefficients)

# Delta-beta (simple group means) for interpretability
ctrl_cols <- which(pheno$case == 0L)
case_cols <- which(pheno$case == 1L)
beta_mean_ctrl <- rowMeans(beta_mat[, ctrl_cols, drop = FALSE])
beta_mean_case <- rowMeans(beta_mat[, case_cols, drop = FALSE])
delta_beta_mean <- beta_mean_case - beta_mean_ctrl

# Model-implied delta-beta via transforming M (using control mean M as baseline)
M_mean_ctrl <- rowMeans(M[, ctrl_cols, drop = FALSE])
coef_case <- fit$coefficients[, coef_idx]
beta_ctrl_pred <- beta_from_M(M_mean_ctrl)
beta_case_pred <- beta_from_M(M_mean_ctrl + coef_case)
delta_beta_model <- beta_case_pred - beta_ctrl_pred

# Delta-method SE on beta scale: d/d(deltaM) logistic2(M0 + deltaM) * SE(deltaM)
se_beta_model <- abs(log(2) * beta_case_pred * (1 - beta_case_pred) * se_M)

res <- data.frame(
  probe_id = rownames(tt),
  logFC_M = tt$logFC,
  se_M = se_M[rownames(tt)],
  t = tt$t,
  P.Value = tt$P.Value,
  adj.P.Val = tt$adj.P.Val,
  delta_beta_mean = delta_beta_mean[rownames(tt)],
  delta_beta_model = delta_beta_model[rownames(tt)],
  se_delta_beta_model = se_beta_model[rownames(tt)],
  mean_beta_ctrl = beta_mean_ctrl[rownames(tt)],
  mean_beta_case = beta_mean_case[rownames(tt)],
  stringsAsFactors = FALSE
)

fwrite(res, file = file.path(opt$out_dir, "limma_dmps.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Persist sample sheet used
samples_out <- pheno[, c("barcode", "gsm", "sex", "age_blood", "braak_stage", "braak_num", "case", cell_types)]
fwrite(samples_out, file = file.path(opt$out_dir, "samples_used.tsv"), sep = "\t", quote = FALSE, na = "NA")

# Session info for reproducibility
sink(file.path(opt$out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("[6/6] Done.")
message("Wrote outputs to: ", opt$out_dir)
