#!/usr/bin/env python3
"""Create exact high-PIP variant-level overlap tables between sex and generation analyses."""

from __future__ import annotations

import argparse
import os
import pandas as pd
import numpy as np
from utils import to_bool, priority_sex, priority_generation


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--meta-loci", required=True)
    p.add_argument("--pip-long", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pip-threshold", type=float, default=0.8)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    meta = pd.read_csv(args.meta_loci, sep=None, engine="python")
    meta.columns = meta.columns.str.strip()
    meta = meta.drop_duplicates("meta_locus_id")
    meta["sex_priority_class"] = meta["sex_source_category"].apply(priority_sex)
    meta["generation_priority_class"] = meta["generation_source_category"].apply(priority_generation)

    usecols = ["locus_id", "analysis_label", "variant_id", "PIP", "inside_core_meta_locus"]
    rows = []
    for chunk in pd.read_csv(args.pip_long, sep="\t", usecols=usecols, chunksize=1_000_000):
        chunk["PIP"] = pd.to_numeric(chunk["PIP"], errors="coerce")
        chunk["inside_core_meta_locus"] = to_bool(chunk["inside_core_meta_locus"])
        sub = chunk[chunk["inside_core_meta_locus"] & (chunk["PIP"] >= args.pip_threshold)].copy()
        if len(sub):
            rows.append(sub)
    high = pd.concat(rows, ignore_index=True).rename(columns={"locus_id": "meta_locus_id"})

    flags = (high.groupby(["meta_locus_id", "variant_id"])
             .agg(
                 highPIP_in_Female=("analysis_label", lambda x: "Female" in set(x)),
                 highPIP_in_Male=("analysis_label", lambda x: "Male" in set(x)),
                 highPIP_in_G1=("analysis_label", lambda x: "G1" in set(x)),
                 highPIP_in_G2=("analysis_label", lambda x: "G2" in set(x)),
                 highPIP_in_G3=("analysis_label", lambda x: "G3" in set(x)),
                 highPIP_in_G4=("analysis_label", lambda x: "G4" in set(x)),
                 analyses_with_highPIP=("analysis_label", lambda x: ";".join(sorted(set(x.dropna().astype(str))))),
                 max_PIP_any=("PIP", "max"),
             )
             .reset_index())
    flags["highPIP_in_sex"] = flags["highPIP_in_Female"] | flags["highPIP_in_Male"]
    flags["highPIP_in_generation"] = flags[["highPIP_in_G1", "highPIP_in_G2", "highPIP_in_G3", "highPIP_in_G4"]].any(axis=1)
    flags = flags.merge(meta[["meta_locus_id", "sex_priority_class", "generation_priority_class"]], on="meta_locus_id", how="left")

    gen_variants = flags[flags["highPIP_in_generation"] &
                         flags["generation_priority_class"].isin(["Generation main", "Generation moderated"])].copy()

    def sex_overlap(row):
        if not row["highPIP_in_sex"]:
            return "Not high-PIP in sex"
        if row["sex_priority_class"] in ["Sex main", "Sex moderated"]:
            return row["sex_priority_class"]
        return "Sex high-PIP but no sex GWAS locus"

    gen_variants["sex_exact_variant_overlap"] = gen_variants.apply(sex_overlap, axis=1)
    panel = pd.crosstab(gen_variants["generation_priority_class"], gen_variants["sex_exact_variant_overlap"])
    panel["Total generation core high-PIP variants"] = panel.sum(axis=1)

    panel.to_csv(os.path.join(args.out_dir, "exact_highPIP_variant_overlap_generation_denominator.tsv"), sep="\t")
    gen_variants.to_csv(os.path.join(args.out_dir, "exact_highPIP_variant_overlap_working_table.tsv"), sep="\t", index=False)
    print(panel)


if __name__ == "__main__":
    main()
