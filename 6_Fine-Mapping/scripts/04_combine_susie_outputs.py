#!/usr/bin/env python3
"""Combine SuSiE PIP/summary files and add core-vs-buffer annotations."""

from __future__ import annotations

import argparse
import glob
import os
import pandas as pd
import numpy as np
from utils import to_bool


DEFAULT_LABELS = {
    "sex_female_meta": "Female",
    "sex_male_meta": "Male",
    "generation1_meta": "G1",
    "generation2_meta": "G2",
    "generation3_meta": "G3",
    "generation4_meta": "G4",
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--result-root", required=True)
    p.add_argument("--metal-mode", required=True)
    p.add_argument("--union-file", required=True, help="Buffered union file with locus_id, START, END, core_START, core_END.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--analyses", nargs="+", default=list(DEFAULT_LABELS.keys()))
    p.add_argument("--labels", nargs="*", default=None,
                   help="Optional labels matching analyses, e.g. Female Male G1 G2 G3 G4.")
    p.add_argument("--pip-threshold", type=float, default=0.8)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.labels:
        if len(args.labels) != len(args.analyses):
            raise ValueError("--labels must have same length as --analyses")
        analysis_labels = dict(zip(args.analyses, args.labels))
    else:
        analysis_labels = {a: DEFAULT_LABELS.get(a, a) for a in args.analyses}

    regions = pd.read_csv(args.union_file, sep="\t")
    if "core_START" not in regions.columns:
        regions["core_START"] = regions["START"]
    if "core_END" not in regions.columns:
        regions["core_END"] = regions["END"]
    regions = regions[["locus_id", "CHR", "START", "END", "core_START", "core_END"]].drop_duplicates()

    all_pip, all_summary = [], []
    for analysis_id, label in analysis_labels.items():
        result_base = os.path.join(args.result_root, analysis_id, args.metal_mode)
        pip_files = sorted(glob.glob(os.path.join(result_base, "susie_results_*", "META", "*.pip.tsv")))
        summary_files = sorted(glob.glob(os.path.join(result_base, "susie_results_*", "META", "*.summary.tsv")))
        print(analysis_id, "pip files:", len(pip_files), "summary files:", len(summary_files))

        for f in pip_files:
            locus_id = os.path.basename(f).replace("_META.pip.tsv", "")
            df = pd.read_csv(f, sep="\t")
            df["locus_id"] = locus_id
            df["analysis_id"] = analysis_id
            df["analysis_label"] = label
            all_pip.append(df)

        for f in summary_files:
            locus_id = os.path.basename(f).replace("_META.summary.tsv", "")
            df = pd.read_csv(f, sep="\t")
            df["locus_id"] = locus_id
            df["analysis_id"] = analysis_id
            df["analysis_label"] = label
            all_summary.append(df)

    pip_long = pd.concat(all_pip, ignore_index=True)
    summary_long = pd.concat(all_summary, ignore_index=True)

    pip_long = pip_long.merge(regions, on="locus_id", how="left", suffixes=("", "_region"))
    pip_long["BP"] = pd.to_numeric(pip_long["BP"], errors="coerce")
    for col in ["START", "END", "core_START", "core_END"]:
        pip_long[col] = pd.to_numeric(pip_long[col], errors="coerce")

    pip_long["inside_core_meta_locus"] = (
        (pip_long["BP"] >= pip_long["core_START"]) & (pip_long["BP"] <= pip_long["core_END"])
    )
    pip_long["inside_buffer_only"] = (
        (pip_long["BP"] >= pip_long["START"]) & (pip_long["BP"] <= pip_long["END"]) &
        (~pip_long["inside_core_meta_locus"])
    )
    pip_long["variant_region"] = "outside"
    pip_long.loc[pip_long["inside_core_meta_locus"], "variant_region"] = "core_meta_locus"
    pip_long.loc[pip_long["inside_buffer_only"], "variant_region"] = "buffer_only"
    pip_long["highPIP_0_8"] = pip_long["PIP"] >= args.pip_threshold

    pip_long_out = os.path.join(args.out_dir, "all_analyses_pip_long_with_core_buffer.tsv")
    summary_long_out = os.path.join(args.out_dir, "all_analyses_susie_summary_long.tsv")
    pip_long.to_csv(pip_long_out, sep="\t", index=False)
    summary_long.to_csv(summary_long_out, sep="\t", index=False)

    # Wide PIP matrix
    key_cols = ["locus_id", "variant_id", "CHR", "BP", "REF", "ALT"]
    pip_wide = pip_long.pivot_table(index=key_cols, columns="analysis_label", values="PIP", aggfunc="max").reset_index()
    pip_wide.columns.name = None
    rename = {"Female": "PIP_Female", "Male": "PIP_Male",
              "G1": "PIP_G1", "G2": "PIP_G2", "G3": "PIP_G3", "G4": "PIP_G4"}
    pip_wide = pip_wide.rename(columns=rename)
    pip_wide_out = os.path.join(args.out_dir, "all_analyses_pip_wide.tsv")
    pip_wide.to_csv(pip_wide_out, sep="\t", index=False)

    print("Saved:", pip_long_out)
    print("Saved:", summary_long_out)
    print("Saved:", pip_wide_out)


if __name__ == "__main__":
    main()
