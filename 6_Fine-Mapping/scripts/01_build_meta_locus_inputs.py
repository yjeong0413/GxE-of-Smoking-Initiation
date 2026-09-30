#!/usr/bin/env python3
"""Build meta-locus coordinate and buffered union files.

Inputs:
  1. Sex loci lifted to hg38.
  2. Generation loci lifted to hg38.
  3. Mapping of original loci to collapsed meta_locus_id.
  4. Optional meta-source summary with sex/generation source categories.

Outputs:
  - meta_loci_for_finemapping_with_hg37_hg38_coordinates.tsv
  - meta_loci_finemapping_job_manifest.tsv
  - meta_loci_union_hg38_for_sumstats.tsv
  - meta_loci_union_hg38_buffered_for_sumstats.tsv
"""

from __future__ import annotations

import argparse
import os
import pandas as pd
import numpy as np
from utils import clean_chr_series, ensure_cols


def build_liftover_table(sex_loci: pd.DataFrame, gen_loci: pd.DataFrame) -> pd.DataFrame:
    sex_cols = ["locus_key", "locus_class", "uniqID", "chr", "start", "end",
                "chr_hg38", "start_hg38", "end_hg38"]
    gen_cols = ["locus_key", "locus_class", "uniqID", "generation_contrast", "chr", "start", "end",
                "chr_hg38", "start_hg38", "end_hg38"]
    ensure_cols(sex_loci, sex_cols, "sex loci")
    ensure_cols(gen_loci, gen_cols, "generation loci")

    sex = sex_loci[sex_cols].drop_duplicates().copy()
    sex["analysis"] = "sex"
    sex = sex.rename(columns={
        "uniqID": "original_uniqID", "chr": "CHR_hg37", "start": "start_hg37", "end": "end_hg37",
        "chr_hg38": "CHR_hg38",
    })

    gen = gen_loci[gen_cols].drop_duplicates().copy()
    gen["analysis"] = "generation"
    gen = gen.rename(columns={
        "uniqID": "original_uniqID", "chr": "CHR_hg37", "start": "start_hg37", "end": "end_hg37",
        "chr_hg38": "CHR_hg38",
    })

    return pd.concat([sex, gen], ignore_index=True, sort=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sex-loci-hg38", required=True)
    p.add_argument("--generation-loci-hg38", required=True)
    p.add_argument("--meta-locus-mapping", required=True,
                   help="TSV with meta_locus_id, analysis, locus_key, locus_class.")
    p.add_argument("--meta-source-summary", default=None,
                   help="Optional table with meta_locus_id, sex_source_category, generation_source_category.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--buffer-bp", type=int, default=400000)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    sex_loci = pd.read_csv(args.sex_loci_hg38, sep="\t")
    gen_loci = pd.read_csv(args.generation_loci_hg38, sep="\t")
    mapping = pd.read_csv(args.meta_locus_mapping, sep="\t")
    mapping.columns = mapping.columns.str.strip()

    liftover = build_liftover_table(sex_loci, gen_loci)

    merged = mapping.merge(
        liftover,
        on=["analysis", "locus_key", "locus_class"],
        how="left",
        suffixes=("", "_liftover"),
    )

    for col in ["start_hg37", "end_hg37", "start_hg38", "end_hg38"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    meta = (
        merged.groupby("meta_locus_id")
        .agg(
            CHR_hg37=("CHR_hg37", "first"),
            start_hg37=("start_hg37", "min"),
            end_hg37=("end_hg37", "max"),
            CHR_hg38=("CHR_hg38", "first"),
            start_hg38=("start_hg38", "min"),
            end_hg38=("end_hg38", "max"),
            n_original_loci=("locus_key", "nunique"),
            n_sex_loci=("analysis", lambda x: (x == "sex").sum()),
            n_generation_loci=("analysis", lambda x: (x == "generation").sum()),
        )
        .reset_index()
    )
    meta["length_hg37_bp"] = meta["end_hg37"] - meta["start_hg37"] + 1
    meta["length_hg38_bp"] = meta["end_hg38"] - meta["start_hg38"] + 1

    if args.meta_source_summary:
        source = pd.read_csv(args.meta_source_summary, sep=None, engine="python")
        source.columns = source.columns.str.strip()
        keep = [c for c in ["meta_locus_id", "sex_source_category", "generation_source_category"] if c in source.columns]
        meta = meta.merge(source[keep].drop_duplicates("meta_locus_id"), on="meta_locus_id", how="left")

    meta_out = os.path.join(args.out_dir, "meta_loci_for_finemapping_with_hg37_hg38_coordinates.tsv")
    meta.to_csv(meta_out, sep="\t", index=False)

    # Job manifest
    strata = ["Female", "Male", "G1", "G2", "G3", "G4"]
    manifest = []
    for _, row in meta.iterrows():
        for s in strata:
            manifest.append({
                "meta_locus_id": row["meta_locus_id"],
                "stratum": s,
                "CHR_hg37": row["CHR_hg37"],
                "start_hg37": row["start_hg37"],
                "end_hg37": row["end_hg37"],
                "CHR_hg38": row["CHR_hg38"],
                "start_hg38": row["start_hg38"],
                "end_hg38": row["end_hg38"],
                "length_hg38_mb": row["length_hg38_bp"] / 1e6,
            })
    pd.DataFrame(manifest).to_csv(
        os.path.join(args.out_dir, "meta_loci_finemapping_job_manifest.tsv"),
        sep="\t", index=False
    )

    # Union file for sumstats extraction
    union = meta[["meta_locus_id", "CHR_hg38", "start_hg38", "end_hg38"]].drop_duplicates().copy()
    union = union.rename(columns={
        "meta_locus_id": "locus_id",
        "CHR_hg38": "CHR",
        "start_hg38": "START",
        "end_hg38": "END",
    })
    union["CHR"] = clean_chr_series(union["CHR"])
    union["START"] = pd.to_numeric(union["START"], errors="raise").astype(int)
    union["END"] = pd.to_numeric(union["END"], errors="raise").astype(int)
    union["locus_size_bp"] = union["END"] - union["START"] + 1
    union = union.sort_values(["CHR", "START", "END"]).reset_index(drop=True)
    union.to_csv(os.path.join(args.out_dir, "meta_loci_union_hg38_for_sumstats.tsv"), sep="\t", index=False)

    buffered = union.copy()
    buffered["core_START"] = buffered["START"]
    buffered["core_END"] = buffered["END"]
    buffered["core_locus_size_bp"] = buffered["core_END"] - buffered["core_START"] + 1
    buffered["START"] = (buffered["core_START"] - args.buffer_bp).clip(lower=1)
    buffered["END"] = buffered["core_END"] + args.buffer_bp
    buffered["locus_size_bp"] = buffered["END"] - buffered["START"] + 1
    buffered["buffer_bp"] = args.buffer_bp
    buffered["buffer_added"] = True
    buffered["core_is_1bp"] = buffered["core_locus_size_bp"] == 1
    buffered["core_is_small_lt_100kb"] = buffered["core_locus_size_bp"] < 100000
    buffered["core_is_small_lt_400kb"] = buffered["core_locus_size_bp"] < 400000
    buffered["buffer_to_core_ratio"] = (
        (buffered["locus_size_bp"] - buffered["core_locus_size_bp"]) / buffered["core_locus_size_bp"]
    )
    buffered["buffer_sensitive_locus"] = buffered["core_is_1bp"] | buffered["core_is_small_lt_100kb"]
    buffered.to_csv(os.path.join(args.out_dir, "meta_loci_union_hg38_buffered_for_sumstats.tsv"),
                    sep="\t", index=False)

    print(f"Saved: {meta_out}")
    print(f"Meta-loci: {meta['meta_locus_id'].nunique()}")


if __name__ == "__main__":
    main()
