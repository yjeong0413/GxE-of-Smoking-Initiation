#!/usr/bin/env python3
"""Check completion of SuSiE runs by comparing joblists with result files."""

from __future__ import annotations

import argparse
import glob
import os
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--meta-base", required=True, help="Directory containing per-analysis joblists.")
    p.add_argument("--result-root", required=True, help="Directory containing per-analysis result subdirectories.")
    p.add_argument("--metal-mode", required=True, help="Result mode subdirectory name.")
    p.add_argument("--analyses", nargs="+", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    joblist_files = {
        "A": "jobs_A_nsnps_le5000_minimal.tsv",
        "B": "jobs_B_nsnps_5001_10000_minimal.tsv",
        "C_split": "jobs_C_split_minimal.tsv",
    }

    rows = []
    for analysis_id in args.analyses:
        expected_loci = set()
        for fname in joblist_files.values():
            joblist = os.path.join(args.meta_base, analysis_id, "joblists", fname)
            if not os.path.exists(joblist):
                print("Missing joblist:", joblist)
                continue
            df = pd.read_csv(joblist, sep="\t")
            expected_loci.update(df["locus_id"].astype(str))

        result_base = os.path.join(args.result_root, analysis_id, args.metal_mode)
        summary_files = glob.glob(os.path.join(result_base, "susie_results_*", "META", "*.summary.tsv"))
        pip_files = glob.glob(os.path.join(result_base, "susie_results_*", "META", "*.pip.tsv"))

        finished_loci = set()
        for f in summary_files:
            finished_loci.add(os.path.basename(f).replace("_META.summary.tsv", ""))

        missing = sorted(expected_loci - finished_loci)
        rows.append({
            "analysis_id": analysis_id,
            "expected_loci": len(expected_loci),
            "summary_files": len(summary_files),
            "pip_files": len(pip_files),
            "finished_loci": len(finished_loci),
            "missing_loci": len(missing),
            "example_missing": ";".join(missing[:5]),
        })

    qc = pd.DataFrame(rows)
    qc.to_csv(args.out, sep="\t", index=False)
    print(qc)
    print("Saved:", args.out)


if __name__ == "__main__":
    main()
