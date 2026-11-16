#!/usr/bin/env python3
"""Preprocess GSE59685 whole-blood methylation data.

This script:
1. Parses GEO metadata to find the 531 study barcodes and the 80 whole-blood barcodes.
2. Filters the signal-intensity file to those samples.
3. Applies sample-level QC (detection failure rate and low-intensity outliers).
4. Applies probe-level QC (detection failure, SNP-affected probes, sex chromosomes,
   Bowtie2 multi-mappers, Chen non-specific probes, and low-variability probes by MAD(M)).
5. Writes sample/probe QC tables and a JSON summary.

The script is designed to be readable and reproducible rather than minimal.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_DATA_DIR = Path("data")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC-preprocess GSE59685 blood samples.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=Path("results/gse59685_qc"))
    parser.add_argument("--detection-threshold", type=float, default=0.01)
    parser.add_argument("--sample-fail-fraction", type=float, default=0.01)
    parser.add_argument("--probe-fail-fraction", type=float, default=0.05)
    parser.add_argument("--variable-fraction", type=float, default=0.10)
    parser.add_argument("--intensity-outlier-iqr", type=float, default=1.5)
    parser.add_argument("--m-offset", type=float, default=100.0)
    parser.add_argument("--chunksize", type=int, default=5000)
    return parser.parse_args()


def parse_sample_metadata(family_path: Path) -> tuple[list[str], list[str], pd.DataFrame]:
    study_rows: list[dict[str, str]] = []
    current: list[str] = []

    with gzip.open(family_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                if current:
                    study_rows.append(_parse_sample_block(current))
                current = [line.rstrip("\n")]
            else:
                current.append(line.rstrip("\n"))

    if current:
        study_rows.append(_parse_sample_block(current))

    sample_df = pd.DataFrame(study_rows)
    study_barcodes = sample_df["barcode"].dropna().astype(str).tolist()
    blood_barcodes = sample_df.loc[sample_df["source_tissue"] == "whole blood", "barcode"].dropna().astype(str).tolist()
    return study_barcodes, blood_barcodes, sample_df


def _parse_sample_block(block: list[str]) -> dict[str, str]:
    sample_id = block[0].split(" = ", 1)[1].strip()
    source = ""
    barcode = ""
    ad_status = ""
    braak = ""
    sex = ""
    age_blood = ""
    age_brain = ""
    title = ""
    for line in block[1:]:
        if line.startswith("!Sample_title = "):
            title = line.split(" = ", 1)[1].strip()
        elif line.startswith("!Sample_source_name_ch1 = "):
            source = line.split(" = ", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = barcode: "):
            barcode = line.split("barcode: ", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = ad.disease.status: "):
            ad_status = line.split(" = ", 1)[1].strip().split(": ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = braak.stage: "):
            braak = line.split(" = ", 1)[1].strip().split(": ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = Sex: "):
            sex = line.split(" = ", 1)[1].strip().split(": ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = age.blood: "):
            age_blood = line.split(" = ", 1)[1].strip().split(": ", 1)[1]
        elif line.startswith("!Sample_characteristics_ch1 = age.brain: "):
            age_brain = line.split(" = ", 1)[1].strip().split(": ", 1)[1]

    return {
        "gsm": sample_id,
        "title": title,
        "source_tissue": source,
        "barcode": barcode,
        "ad_status": ad_status,
        "braak_stage": braak,
        "sex": sex,
        "age_blood": age_blood,
        "age_brain": age_brain,
    }


def read_manifest(manifest_path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(
        manifest_path,
        skiprows=7,
        usecols=["IlmnID", "CHR", "Probe_SNPs_10"],
        dtype=str,
        keep_default_na=False,
    )
    manifest["CHR"] = manifest["CHR"].str.strip()
    manifest["Probe_SNPs_10"] = manifest["Probe_SNPs_10"].str.strip()
    return manifest


def read_probe_list(path: Path) -> set[str]:
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "TargetID" in df.columns:
            values = df["TargetID"].astype(str)
        else:
            values = df.iloc[:, 0].astype(str)
    else:
        values = pd.read_csv(path, header=None, dtype=str, keep_default_na=False).iloc[:, 0].astype(str)
    return {value.strip() for value in values if value and value.strip()}


def signal_header(signal_path: Path) -> list[str]:
    with gzip.open(signal_path, "rt", encoding="utf-8", errors="replace") as handle:
        handle.readline()
        handle.readline()
        header = handle.readline().rstrip("\n").split(",")
    return header


def compute_sample_intensity_qc(
    signal_path: Path,
    blood_barcodes: list[str],
    chunksize: int,
    iqr_multiplier: float,
) -> tuple[pd.DataFrame, list[str]]:
    sample_totals: dict[str, list[np.ndarray]] = {barcode: [] for barcode in blood_barcodes}
    probe_count = 0

    usecols = ["Unnamed: 0"]
    for barcode in blood_barcodes:
        usecols.extend([f"{barcode} Unmethylated Signal", f"{barcode} Methylated Signal"])

    for chunk in pd.read_csv(signal_path, compression="gzip", skiprows=2, usecols=usecols, chunksize=chunksize):
        unmethylated = chunk.iloc[:, 1::2].to_numpy(dtype=float)
        methylated = chunk.iloc[:, 2::2].to_numpy(dtype=float)
        totals = unmethylated + methylated
        for index, barcode in enumerate(blood_barcodes):
            sample_totals[barcode].append(totals[:, index].astype(np.float32, copy=False))
        probe_count += len(chunk)

    rows = []
    for barcode in blood_barcodes:
        values = np.concatenate(sample_totals[barcode]) if sample_totals[barcode] else np.array([], dtype=np.float32)
        med_total = float(np.nanmedian(values))
        mean_total = float(np.nanmean(values))
        rows.append({"barcode": barcode, "median_total_intensity": med_total, "mean_total_intensity": mean_total})

    sample_df = pd.DataFrame(rows)
    q1 = float(sample_df["median_total_intensity"].quantile(0.25))
    q3 = float(sample_df["median_total_intensity"].quantile(0.75))
    iqr = q3 - q1
    threshold = q1 - iqr_multiplier * iqr
    sample_df["low_intensity_outlier"] = sample_df["median_total_intensity"] < threshold
    sample_df["low_intensity_threshold"] = threshold
    sample_df["probe_count"] = probe_count
    return sample_df, [barcode for barcode in sample_df.loc[~sample_df["low_intensity_outlier"], "barcode"].tolist()]


def compute_probe_qc(
    signal_path: Path,
    kept_barcodes: list[str],
    detection_threshold: float,
    chunksize: int,
    m_offset: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    usecols = ["Unnamed: 0"]
    for barcode in kept_barcodes:
        usecols.extend(
            [
                f"{barcode} Unmethylated Signal",
                f"{barcode} Methylated Signal",
                f"{barcode} Detection Pval",
            ]
        )

    probe_ids: list[str] = []
    detection_fail_counts: list[int] = []
    mad_m_values: list[float] = []
    total_probe_count = 0
    sample_fail_counts = np.zeros(len(kept_barcodes), dtype=np.int64)

    for chunk in pd.read_csv(signal_path, compression="gzip", skiprows=2, usecols=usecols, chunksize=chunksize):
        probe_ids.extend(chunk.iloc[:, 0].astype(str).tolist())

        unmethylated = chunk.iloc[:, 1::3].to_numpy(dtype=float)
        methylated = chunk.iloc[:, 2::3].to_numpy(dtype=float)
        detection = chunk.iloc[:, 3::3].to_numpy(dtype=float)

        fail_mask = np.isnan(detection) | (detection > detection_threshold)
        detection_fail_counts.extend(fail_mask.sum(axis=1).astype(int).tolist())
        sample_fail_counts += fail_mask.sum(axis=0).astype(np.int64)

        m_values = np.log2((methylated + m_offset) / (unmethylated + m_offset))
        median_m = np.nanmedian(m_values, axis=1, keepdims=True)
        mad_m = np.nanmedian(np.abs(m_values - median_m), axis=1)
        mad_m_values.extend(mad_m.astype(float).tolist())
        total_probe_count += len(chunk)

    probe_df = pd.DataFrame(
        {
            "probe_id": probe_ids,
            "detection_fail_count": detection_fail_counts,
            "detection_fail_fraction": np.array(detection_fail_counts, dtype=float) / len(kept_barcodes),
            "mad_m": mad_m_values,
        }
    )

    sample_df = pd.DataFrame(
        {
            "barcode": kept_barcodes,
            "detection_fail_count": sample_fail_counts,
            "detection_fail_fraction": sample_fail_counts / total_probe_count,
            "probe_count": total_probe_count,
        }
    )
    sample_df["high_detection_fail_sample"] = sample_df["detection_fail_fraction"] > 0.01
    return probe_df, sample_df


def annotate_probes(
    probe_df: pd.DataFrame,
    manifest: pd.DataFrame,
    chen_blacklist: set[str],
    bowtie_blacklist: set[str],
    variable_fraction: float,
) -> pd.DataFrame:
    manifest_index = manifest.set_index("IlmnID")
    annotated = probe_df.copy()
    annotated["in_manifest"] = annotated["probe_id"].isin(manifest_index.index)

    manifest_rows = manifest_index.reindex(annotated["probe_id"])
    annotated["chr"] = manifest_rows["CHR"].to_numpy()
    annotated["probe_snp_10"] = manifest_rows["Probe_SNPs_10"].to_numpy()

    chr_upper = annotated["chr"].fillna("").str.upper()
    annotated["sex_chr"] = chr_upper.isin({"X", "Y", "CHRX", "CHRY"})
    annotated["snp10"] = annotated["probe_snp_10"].fillna("").astype(str).str.strip().ne("")
    annotated["bowtie_multimap"] = annotated["probe_id"].isin(bowtie_blacklist)
    annotated["chen_non_specific"] = annotated["probe_id"].isin(chen_blacklist)
    annotated["cross_reactive_union"] = annotated["bowtie_multimap"] | annotated["chen_non_specific"]
    annotated["technical_keep"] = ~(
        annotated["snp10"] | annotated["sex_chr"] | annotated["cross_reactive_union"] | annotated["detection_fail_fraction"].gt(0.05)
    )

    technical_pass = annotated.loc[annotated["technical_keep"], "mad_m"]
    variable_threshold = float(technical_pass.quantile(variable_fraction))
    annotated["low_variability"] = annotated["mad_m"] <= variable_threshold
    annotated["final_keep"] = annotated["technical_keep"] & ~annotated["low_variability"]
    annotated["mad_threshold"] = variable_threshold
    return annotated


def write_outputs(
    output_dir: Path,
    sample_df: pd.DataFrame,
    probe_df: pd.DataFrame,
    summary: dict,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_df.to_csv(output_dir / "sample_qc.tsv.gz", sep="\t", index=False, compression="gzip")
    probe_df.to_csv(output_dir / "probe_qc.tsv.gz", sep="\t", index=False, compression="gzip")
    with open(output_dir / "qc_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


def main() -> None:
    args = parse_args()

    
    family_path = args.data_dir / "GSE59685_family.soft.gz"
    signal_path = args.data_dir / "GSE59685_signal_intensities.csv.gz"
    manifest_path = args.data_dir / "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
    chen_path = args.data_dir / "chen_48639_non_specific_450k.csv"
    bowtie_path = args.data_dir / "bowtie_multimap_450k.txt"

    study_barcodes, blood_barcodes, sample_meta = parse_sample_metadata(family_path)
    header = signal_header(signal_path)
    study_barcodes = [barcode for barcode in study_barcodes if f"{barcode} Unmethylated Signal" in header]
    blood_barcodes = [barcode for barcode in blood_barcodes if f"{barcode} Unmethylated Signal" in header]

    manifest = read_manifest(manifest_path)
    chen_blacklist = read_probe_list(chen_path)
    bowtie_blacklist = read_probe_list(bowtie_path)

    sample_intensity_df, kept_barcodes = compute_sample_intensity_qc(
        signal_path=signal_path,
        blood_barcodes=blood_barcodes,
        chunksize=args.chunksize,
        iqr_multiplier=args.intensity_outlier_iqr,
    )

    probe_df, sample_detection_df = compute_probe_qc(
        signal_path=signal_path,
        kept_barcodes=kept_barcodes,
        detection_threshold=args.detection_threshold,
        chunksize=args.chunksize,
        m_offset=args.m_offset,
    )

    annotated_probe_df = annotate_probes(
        probe_df=probe_df,
        manifest=manifest,
        chen_blacklist=chen_blacklist,
        bowtie_blacklist=bowtie_blacklist,
        variable_fraction=args.variable_fraction,
    )

    sample_df = sample_intensity_df.merge(sample_detection_df, on="barcode", how="left", suffixes=("", "_detection"))
    sample_df = sample_df.merge(sample_meta, on="barcode", how="left")
    sample_df["sample_detection_fail_keep"] = sample_df["detection_fail_fraction"] <= args.sample_fail_fraction
    sample_df["sample_final_keep"] = sample_df["sample_detection_fail_keep"] & ~sample_df["low_intensity_outlier"]

    final_probe_count = int(annotated_probe_df["final_keep"].sum())
    summary = {
        "study_barcodes_found": len(study_barcodes),
        "blood_barcodes_found": len(blood_barcodes),
        "sample_low_intensity_outliers": int(sample_df["low_intensity_outlier"].sum()),
        "sample_high_detection_fail_count": int(sample_df["sample_detection_fail_keep"].eq(False).sum()),
        "probe_total": int(len(annotated_probe_df)),
        "probe_detection_fail_gt_5pct": int((annotated_probe_df["detection_fail_fraction"] > args.probe_fail_fraction).sum()),
        "probe_snp10_count": int(annotated_probe_df["snp10"].sum()),
        "probe_sex_chr_count": int(annotated_probe_df["sex_chr"].sum()),
        "probe_bowtie_multimap_count": int(annotated_probe_df["bowtie_multimap"].sum()),
        "probe_chen_non_specific_count": int(annotated_probe_df["chen_non_specific"].sum()),
        "probe_cross_reactive_union_count": int(annotated_probe_df["cross_reactive_union"].sum()),
        "probe_low_variability_count": int(annotated_probe_df["low_variability"].sum()),
        "probe_final_keep_count": final_probe_count,
        "chen_blacklist_size": len(chen_blacklist),
        "bowtie_blacklist_size": len(bowtie_blacklist),
        "blacklist_overlap_count": len(chen_blacklist & bowtie_blacklist),
        "blacklist_union_count": len(chen_blacklist | bowtie_blacklist),
        "blacklist_union_overlap_with_matrix": int(annotated_probe_df["cross_reactive_union"].sum()),
        "mad_threshold": float(annotated_probe_df["mad_threshold"].iloc[0]),
        "detection_threshold": args.detection_threshold,
        "sample_fail_fraction": args.sample_fail_fraction,
        "probe_fail_fraction": args.probe_fail_fraction,
        "variable_fraction": args.variable_fraction,
        "intensity_outlier_iqr": args.intensity_outlier_iqr,
        "m_offset": args.m_offset,
        "sample_output_barcodes": len(kept_barcodes),
    }

    write_outputs(args.output_dir, sample_df, annotated_probe_df, summary)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
