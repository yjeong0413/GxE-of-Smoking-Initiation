#!/usr/bin/env python3
"""Create A/B/C joblists from a locus sumstats extraction summary."""

from __future__ import annotations

import argparse
import os
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, help="locus_sumstats_extraction_summary.tsv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--a-max-snps", type=int, default=5000)
    p.add_argument("--b-max-snps", type=int, default=10000)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    summary = pd.read_csv(args.summary, sep="\t")
    summary = summary.loc[summary["status"] == "kept"].copy()

    a = summary.loc[summary["n_snps"] <= args.a_max_snps].copy()
    a[["ancestry", "locus_id", "CHR", "n_snps", "sumstats_file", "snplist_file"]].to_csv(
        os.path.join(args.out_dir, "jobs_A_nsnps_le5000_minimal.tsv"), sep="\t", index=False
    )

    b = summary.loc[(summary["n_snps"] > args.a_max_snps) & (summary["n_snps"] <= args.b_max_snps)].copy()
    b[["ancestry", "locus_id", "CHR", "n_snps", "sumstats_file", "snplist_file"]].to_csv(
        os.path.join(args.out_dir, "jobs_B_nsnps_5001_10000_minimal.tsv"), sep="\t", index=False
    )

    c = summary.loc[summary["n_snps"] > args.b_max_snps].copy()
    c["original_locus_id"] = c["locus_id"]
    c[["ancestry", "original_locus_id", "locus_id", "CHR", "START", "END",
       "n_snps", "sumstats_file", "snplist_file"]].to_csv(
        os.path.join(args.out_dir, "jobs_C_split_minimal.tsv"), sep="\t", index=False
    )

    print(f"A jobs: {len(a)}")
    print(f"B jobs: {len(b)}")
    print(f"C jobs: {len(c)}")


if __name__ == "__main__":
    main()
