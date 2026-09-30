#!/usr/bin/env python3
"""Create source, fine-mapping, and UpSet membership tables for meta-loci."""

from __future__ import annotations

import argparse
import os
import pandas as pd
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--meta-loci", required=True, help="Meta-locus table with source categories.")
    p.add_argument("--core-buffer-summary", required=True, help="meta_locus_core_buffer_highPIP_summary_by_analysis.tsv")
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    meta = pd.read_csv(args.meta_loci, sep=None, engine="python")
    meta.columns = meta.columns.str.strip()

    meta["has_sex_main"] = meta["sex_source_category"].isin(["sex_main only", "sex_main + sex_moderated"])
    meta["has_sex_moderated"] = meta["sex_source_category"].isin(["sex_moderated only", "sex_main + sex_moderated"])
    meta["has_generation_main"] = meta["generation_source_category"].isin(["generation_main only", "generation_main + generation_moderated"])
    meta["has_generation_moderated"] = meta["generation_source_category"].isin(["generation_moderated only", "generation_main + generation_moderated"])

    core = pd.read_csv(args.core_buffer_summary, sep="\t")
    core.columns = core.columns.str.strip()
    core["has_highPIP_core"] = core["has_highPIP_core"].astype(bool)

    core_wide = core.pivot_table(index="locus_id", columns="analysis_label",
                                 values="has_highPIP_core", aggfunc="max", fill_value=False).reset_index()
    core_wide = core_wide.rename(columns={
        "locus_id": "meta_locus_id",
        "Female": "FM_core_Female",
        "Male": "FM_core_Male",
        "G1": "FM_core_G1",
        "G2": "FM_core_G2",
        "G3": "FM_core_G3",
        "G4": "FM_core_G4",
    })
    for c in ["FM_core_Female", "FM_core_Male", "FM_core_G1", "FM_core_G2", "FM_core_G3", "FM_core_G4"]:
        if c not in core_wide.columns:
            core_wide[c] = False
        core_wide[c] = core_wide[c].astype(bool)

    core_wide["has_sex_finemap_core"] = core_wide["FM_core_Female"] | core_wide["FM_core_Male"]
    core_wide["has_generation_finemap_core"] = (
        core_wide["FM_core_G1"] | core_wide["FM_core_G2"] | core_wide["FM_core_G3"] | core_wide["FM_core_G4"]
    )
    core_wide["n_generation_core_highPIP"] = core_wide[["FM_core_G1", "FM_core_G2", "FM_core_G3", "FM_core_G4"]].sum(axis=1)

    meta = meta.merge(core_wide, on="meta_locus_id", how="left")
    bool_cols = [c for c in meta.columns if c.startswith("FM_core_")] + ["has_sex_finemap_core", "has_generation_finemap_core"]
    for c in bool_cols:
        meta[c] = meta[c].fillna(False).astype(bool)
    meta["n_generation_core_highPIP"] = meta["n_generation_core_highPIP"].fillna(0).astype(int)

    meta_out = os.path.join(args.out_dir, "meta_locus_source_and_core_finemap_evidence.tsv")
    meta.to_csv(meta_out, sep="\t", index=False)

    rows = []
    defs = [
        ("sex_main", "sex", "has_sex_main", "has_sex_finemap_core"),
        ("sex_moderated", "sex", "has_sex_moderated", "has_sex_finemap_core"),
        ("generation_main", "generation", "has_generation_main", "has_generation_finemap_core"),
        ("generation_moderated", "generation", "has_generation_moderated", "has_generation_finemap_core"),
    ]
    for effect_class, domain, source_col, fm_col in defs:
        sub = meta[meta[source_col]]
        rows.append({
            "effect_class": effect_class,
            "domain": domain,
            "n_meta_loci": sub["meta_locus_id"].nunique(),
            "n_core_finemapped_meta_loci": int(sub[fm_col].sum()),
            "n_not_core_finemapped_meta_loci": int((~sub[fm_col]).sum()),
            "percent_core_finemapped": round(100 * sub[fm_col].mean(), 2) if len(sub) else np.nan,
        })
    pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, "main_moderated_effect_level_core_finemap_summary.tsv"),
                              sep="\t", index=False)

    source_matrix = pd.crosstab(meta["generation_source_category"], meta["sex_source_category"])
    source_matrix.to_csv(os.path.join(args.out_dir, "sex_generation_source_overlap_matrix.tsv"), sep="\t")

    upset = meta[[
        "meta_locus_id", "has_sex_main", "has_sex_moderated",
        "has_generation_main", "has_generation_moderated",
        "has_sex_finemap_core", "has_generation_finemap_core",
        "FM_core_Female", "FM_core_Male", "FM_core_G1", "FM_core_G2", "FM_core_G3", "FM_core_G4",
    ]].copy().rename(columns={
        "has_sex_main": "Sex_main_GWAS",
        "has_sex_moderated": "Sex_moderated_GWAS",
        "has_generation_main": "Generation_main_GWAS",
        "has_generation_moderated": "Generation_moderated_GWAS",
        "has_sex_finemap_core": "Sex_core_finemap",
        "has_generation_finemap_core": "Generation_core_finemap",
    })
    upset.to_csv(os.path.join(args.out_dir, "upset_membership_table_core_finemap.tsv"), sep="\t", index=False)

    set_cols = ["Sex_main_GWAS", "Sex_moderated_GWAS", "Generation_main_GWAS",
                "Generation_moderated_GWAS", "Sex_core_finemap", "Generation_core_finemap"]
    combo = upset.groupby(set_cols).size().reset_index(name="n_meta_loci").sort_values("n_meta_loci", ascending=False)
    combo["combination"] = combo.apply(lambda r: " + ".join([c for c in set_cols if bool(r[c])]) or "None", axis=1)
    combo[["combination", "n_meta_loci"] + set_cols].to_csv(
        os.path.join(args.out_dir, "upset_combination_counts_core_finemap.tsv"), sep="\t", index=False
    )

    print("Saved:", meta_out)


if __name__ == "__main__":
    main()
