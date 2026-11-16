# GSE59685 Whole-Blood 450K Methylation Analysis

This repository contains a reproducible pipeline for the take-home task:

- Dataset: GEO `GSE59685` (Illumina HumanMethylation450 BeadChip; 531 arrays total).
- Analysis subset: whole blood only (80 arrays; 72 used in the final model after Braak binning).
- Goal: differential methylation (cases vs controls) after converting Braak stage to a binary variable.

The repo tracks code, QC outputs, and plots. Large raw inputs live in `data/` and are intentionally **not** tracked in git.

## Final Result

The main end-to-end output of this repository is the written report: `reports/FinalReport.pdf`.

## Repo Layout

- `data_summary.txt`: a short, plain-English walkthrough of what the GEO files contain and what each pipeline step does.
- `scripts/01_preprocess_gse59685.py`: blood-only QC + probe/sample filtering.
- `scripts/02_differential_analysis.R`: Houseman-style cell composition adjustment + limma DMP testing on M-values.
- `scripts/03_plot_results.R`: annotation + tables + plots (volcano/QQ/p-hist/boxplots/heatmap/manhattan/regional).
- `results/`: QC tables, limma results, and plot outputs.
- `reports/`: LaTeX report source and PDFs.

## Quickstart

1. Download data into `data/` (not tracked):
   ```bash
   bash scripts/00_download_geo_data.sh
   ```

2. Run QC / filtering:
   ```bash
   python3 scripts/01_preprocess_gse59685.py
   ```

3. Run differential methylation:
   ```bash
   Rscript scripts/02_differential_analysis.R
   ```

4. Generate plots + annotated tables:
   ```bash
   Rscript scripts/03_plot_results.R
   ```

Outputs:
- QC summaries: `results/gse59685_qc/`
- Differential results: `results/differential_analysis/`
- Plots/tables: `results/differential_analysis/plots/`

## Installation

### Option A (Recommended): Conda/Mamba Environment

This avoids the common macOS issues compiling Bioconductor packages like `minfi`.

```bash
# Install micromamba/mamba/conda however you prefer, then:
mamba create -n gse59685 \
  -c conda-forge -c bioconda \
  python=3.11 python-numpy python-pandas \
  r-base r-optparse r-data.table r-quadprog \
  bioconductor-limma bioconductor-minfi bioconductor-flowsorted.blood.450k bioconductor-summarizedexperiment

mamba activate gse59685
```

### Option B: System Python + System R

Python:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install numpy pandas
```

R (CRAN + Bioconductor):
```r
install.packages("BiocManager", repos="https://cloud.r-project.org")
options(repos = BiocManager::repositories())
BiocManager::install(c("limma","minfi","FlowSorted.Blood.450k","SummarizedExperiment"), ask=FALSE, update=TRUE)
install.packages(c("optparse","data.table","quadprog"), repos="https://cloud.r-project.org")
```

## Data: How We Handle It In Git

### Why `data/` is not committed

The key GEO inputs are huge (e.g. `GSE59685_betas.csv.gz` is ~1.9 GB; `signal_intensities` is ~1.5 GB). GitHub rejects regular git pushes with files >100 MB, so pushing raw data directly to git is not practical.

Instead:
- `data/` is ignored via `.gitignore`.
- `scripts/00_download_geo_data.sh` re-downloads the exact GEO files and the filtering blacklists.

### Can we push the data directly to git?

Not with normal git on GitHub (the files are far above the 100 MB per-file hard limit).

If you absolutely need “data tracked alongside code”, options are:
- **Git LFS**: works technically, but large files quickly burn through LFS storage/bandwidth quotas (and requires contributors to install `git-lfs`).
- **External object storage** (recommended for big data): S3 / GCS / institutional storage, with a small `download.sh` in the repo.
- **Data versioning tools**: DVC or DataLad/git-annex (store pointers in git; data in remote storage).
- **GitHub Releases**: can host assets, but still has size limits and isn’t ideal for multi-GB datasets.

For a public GEO dataset like this, the simplest approach is to keep raw inputs out of git and download on demand.

## Notes

- `results/differential_analysis/plots/annotated_results.tsv` is intentionally ignored (it is >100 MB and will be rejected by GitHub). Regenerate it with `Rscript scripts/03_plot_results.R`.
- `results/differential_analysis/limma_dmps.tsv` is ~69 MB; GitHub allows it, but may warn about recommended size limits.
