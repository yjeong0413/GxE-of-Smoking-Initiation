#!/usr/bin/env python3
"""Summarize core high-PIP patterns and buffer-only sensitivity."""

from __future__ import annotations

import argparse
import os
import pandas as pd
import numpy as np
from utils import to_bool


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pip-wide", required=True)
    p.add_argument("--pip-long", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pip-threshold", type=float, default=0.8)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    pip = pd.read_csv(args.pip_wide, sep="\t")
    for c in ["PIP_Female", "PIP_Male", "PIP_G1", "PIP_G2", "PIP_G3", "PIP_G4"]:
        if c not in pip.columns:
            pip[c] = np.nan

    summary = (
        pip.groupby("locus_id")
        .agg(
            n_variants=("variant_id", "nunique"),
            max_PIP_Female=("PIP_Female", "max"),
            max_PIP_Male=("PIP_Male", "max"),
            max_PIP_G1=("PIP_G1", "max"),
            max_PIP_G2=("PIP_G2", "max"),
            max_PIP_G3=("PIP_G3", "max"),
            max_PIP_G4=("PIP_G4", "max"),
        )
        .reset_index()
    )
    t = args.pip_threshold
    summary["sex_abs_PIP_difference"] = (summary["max_PIP_Female"] - summary["max_PIP_Male"]).abs()
    summary["sex_highPIP_pattern"] = np.select(
        [
            (summary["max_PIP_Female"] >= t) & (summary["max_PIP_Male"] >= t),
            (summary["max_PIP_Female"] >= t) & (summary["max_PIP_Male"] < t),
            (summary["max_PIP_Female"] < t) & (summary["max_PIP_Male"] >= t),
        ],
        ["Shared high-PIP", "Female high-PIP only", "Male high-PIP only"],
        default="No high-PIP",
    )
    gen_cols = ["max_PIP_G1", "max_PIP_G2", "max_PIP_G3", "max_PIP_G4"]
    summary["max_PIP_any_generation"] = summary[gen_cols].max(axis=1)
    summary["min_PIP_any_generation"] = summary[gen_cols].min(axis=1)
    summary["generation_PIP_range"] = summary["max_PIP_any_generation"] - summary["min_PIP_any_generation"]
    summary["n_generations_highPIP"] = (summary[gen_cols] >= t).sum(axis=1)
    summary["generation_highPIP_pattern"] = np.select(
        [summary["n_generations_highPIP"] >= 2, summary["n_generations_highPIP"] == 1],
        ["Shared across generations", "Generation high-PIP only"],
        default="No high-PIP",
    )
    out = os.path.join(args.out_dir, "meta_locus_highPIP_summary.tsv")
    summary.to_csv(out, sep="\t", index=False)

    long = pd.read_csv(args.pip_long, sep="\t")
    long["PIP"] = pd.to_numeric(long["PIP"], errors="coerce")
    long["inside_core_meta_locus"] = to_bool(long["inside_core_meta_locus"])
    long["inside_buffer_only"] = to_bool(long["inside_buffer_only"])
    long["highPIP"] = long["PIP"] >= t

    core = (long[long["inside_core_meta_locus"]]
            .groupby(["analysis_id", "analysis_label", "locus_id"])
            .agg(max_PIP_core=("PIP", "max"), n_highPIP_core=("highPIP", "sum"))
            .reset_index())
    buffer = (long[long["inside_buffer_only"]]
              .groupby(["analysis_id", "analysis_label", "locus_id"])
              .agg(max_PIP_buffer_only=("PIP", "max"), n_highPIP_buffer_only=("highPIP", "sum"))
              .reset_index())
    anywhere = (long.groupby(["analysis_id", "analysis_label", "locus_id"])
                .agg(n_variants=("variant_id", "nunique"), max_PIP_anywhere=("PIP", "max"))
                .reset_index())
    qc = anywhere.merge(core, how="left").merge(buffer, how="left")
    for c in ["max_PIP_core", "max_PIP_buffer_only"]:
        qc[c] = qc[c].fillna(0)
    for c in ["n_highPIP_core", "n_highPIP_buffer_only"]:
        qc[c] = qc[c].fillna(0).astype(int)
    qc["has_highPIP_anywhere"] = qc["max_PIP_anywhere"] >= t
    qc["has_highPIP_core"] = qc["max_PIP_core"] >= t
    qc["has_highPIP_buffer_only"] = qc["max_PIP_buffer_only"] >= t
    qc["highPIP_region_class"] = np.select(
        [
            qc["has_highPIP_core"] & qc["has_highPIP_buffer_only"],
            qc["has_highPIP_core"] & ~qc["has_highPIP_buffer_only"],
            ~qc["has_highPIP_core"] & qc["has_highPIP_buffer_only"],
        ],
        ["core_and_buffer_highPIP", "core_only_highPIP", "buffer_only_highPIP"],
        default="no_highPIP",
    )
    qc_out = os.path.join(args.out_dir, "meta_locus_core_buffer_highPIP_summary_by_analysis.tsv")
    qc.to_csv(qc_out, sep="\t", index=False)
    print("Saved:", out)
    print("Saved:", qc_out)


if __name__ == "__main__":
    main()
