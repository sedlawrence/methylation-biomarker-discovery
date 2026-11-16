#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="${DATA_DIR:-data}"
mkdir -p "$DATA_DIR"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "[skip] $out"
    return 0
  fi
  echo "[get]  $url"
  curl -L --fail --retry 3 --connect-timeout 30 -o "$out" "$url"
}

GSE="GSE59685"
GSE_SERIES_DIR="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/${GSE}"
SUPPL="${GSE_SERIES_DIR}/suppl"
SOFT="${GSE_SERIES_DIR}/soft"
MATRIX="${GSE_SERIES_DIR}/matrix"

# Main study files
fetch "${SUPPL}/${GSE}_betas.csv.gz" "${DATA_DIR}/${GSE}_betas.csv.gz"
fetch "${SUPPL}/${GSE}_signal_intensities.csv.gz" "${DATA_DIR}/${GSE}_signal_intensities.csv.gz"
fetch "${SUPPL}/${GSE}_SampleBarcode_Re-analyzedGSMs_NewGSMs.txt.gz" "${DATA_DIR}/${GSE}_SampleBarcode_Re-analyzedGSMs_NewGSMs.txt.gz"

# Metadata (used by scripts/01_preprocess_gse59685.py)
fetch "${SOFT}/${GSE}_family.soft.gz" "${DATA_DIR}/${GSE}_family.soft.gz"
fetch "${MATRIX}/${GSE}_series_matrix.txt.gz" "${DATA_DIR}/${GSE}_series_matrix.txt.gz"

# Platform/annotation bundle distributed as supplementary "RAW"
fetch "${SUPPL}/${GSE}_RAW.tar" "${DATA_DIR}/${GSE}_RAW.tar"
if [[ ! -s "${DATA_DIR}/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz" ]]; then
  echo "[tar]  extracting platform manifest files from ${DATA_DIR}/${GSE}_RAW.tar"
  tar -xf "${DATA_DIR}/${GSE}_RAW.tar" -C "${DATA_DIR}"
fi

# Cross-reactive/multi-mapping probe lists (sirselim/illumina450k_filtering)
FILTERING_BASE="https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master"
fetch "${FILTERING_BASE}/48639-non-specific-probes-Illumina450k.csv" "${DATA_DIR}/chen_48639_non_specific_450k.csv"
fetch "${FILTERING_BASE}/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt" "${DATA_DIR}/bowtie_multimap_450k.txt"
fetch "${FILTERING_BASE}/README.md" "${DATA_DIR}/illumina450k_filtering_README.md"

echo "[done] data downloaded to: ${DATA_DIR}/"

