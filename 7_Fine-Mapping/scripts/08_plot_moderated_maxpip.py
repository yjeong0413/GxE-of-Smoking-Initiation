#!/usr/bin/env python3
"""Plot max-PIP patterns across sex- and generation-moderated meta-loci."""

from __future__ import annotations

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import to_bool


def format_mb(x):
    if pd.isna(x):
        return "NA"
    return f"{float(x) / 1e6:.1f}"


def make_locus_label(df):
    labels = []
    for _, r in df.iterrows():
        chrom = str(r.get("chr_hg38", r.get("CHR_hg38", "NA"))).replace("chr", "")
        start = format_mb(r.get("meta_start_hg38", r.get("start_hg38", np.nan)))
        end = format_mb(r.get("meta_end_hg38", r.get("end_hg38", np.nan)))
        labels.append(f"{r['meta_locus_id']} | chr{chrom}:{start}-{end}Mb")
    return labels


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--moderated-summary", required=True,
                   help="Table with meta_locus_id, moderated_domain, coordinates.")
    p.add_argument("--pip-long", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pip-threshold", type=float, default=0.8)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    summary = pd.read_csv(args.moderated_summary, sep="\t")
    loci = set(summary["meta_locus_id"].astype(str))

    rows = []
    usecols = ["locus_id", "analysis_label", "PIP", "inside_core_meta_locus"]
    for chunk in pd.read_csv(args.pip_long, sep="\t", usecols=usecols, chunksize=1_000_000):
        chunk["locus_id"] = chunk["locus_id"].astype(str)
        chunk["PIP"] = pd.to_numeric(chunk["PIP"], errors="coerce")
        chunk["inside_core_meta_locus"] = to_bool(chunk["inside_core_meta_locus"])
        sub = chunk[chunk["locus_id"].isin(loci) & chunk["inside_core_meta_locus"]].copy()
        if len(sub):
            rows.append(sub)
    pip_core = pd.concat(rows, ignore_index=True)
    max_pip = pip_core.groupby(["locus_id", "analysis_label"])["PIP"].max().reset_index().rename(columns={"locus_id": "meta_locus_id"})

    # Sex plot
    sex_loci = summary[summary["moderated_domain"].eq("sex_moderated")].copy()
    sex_wide = max_pip[max_pip["analysis_label"].isin(["Female", "Male"])].pivot_table(
        index="meta_locus_id", columns="analysis_label", values="PIP", aggfunc="max"
    ).reset_index().rename(columns={"Female": "PIP_Female", "Male": "PIP_Male"})
    sex_plot = sex_loci.merge(sex_wide, on="meta_locus_id", how="left")
    for c in ["PIP_Female", "PIP_Male"]:
        sex_plot[c] = pd.to_numeric(sex_plot.get(c, np.nan), errors="coerce").fillna(0)
    sex_plot["max_any_PIP"] = sex_plot[["PIP_Female", "PIP_Male"]].max(axis=1)
    sex_plot["locus_label"] = make_locus_label(sex_plot)
    sex_plot = sex_plot.sort_values(["max_any_PIP", "meta_locus_id"]).reset_index(drop=True)

    plt.figure(figsize=(8.5, max(8, 0.22 * len(sex_plot))))
    y = np.arange(len(sex_plot))
    for i, r in sex_plot.iterrows():
        plt.plot([r["PIP_Female"], r["PIP_Male"]], [i, i], linewidth=1, alpha=0.6)
    plt.scatter(sex_plot["PIP_Female"], y, label="Female", s=28)
    plt.scatter(sex_plot["PIP_Male"], y, label="Male", s=28)
    plt.axvline(args.pip_threshold, linestyle=":", linewidth=1)
    plt.yticks(y, sex_plot["locus_label"], fontsize=7)
    plt.xlabel("Maximum PIP within core meta-locus")
    plt.ylabel("Sex-moderated meta-loci")
    plt.title("Female and male max-PIP patterns across\nsex-moderated meta-loci")
    plt.xlim(-0.02, 1.03)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "sex_moderated_maxPIP_core_plot.png"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "sex_moderated_maxPIP_core_plot.pdf"))
    plt.close()
    sex_plot.to_csv(os.path.join(args.out_dir, "sex_moderated_maxPIP_core_plot_data.tsv"), sep="\t", index=False)

    # Generation plot
    gen_loci = summary[summary["moderated_domain"].eq("generation_moderated")].copy()
    gen_wide = max_pip[max_pip["analysis_label"].isin(["G1", "G2", "G3", "G4"])].pivot_table(
        index="meta_locus_id", columns="analysis_label", values="PIP", aggfunc="max"
    ).reset_index().rename(columns={"G1": "PIP_G1", "G2": "PIP_G2", "G3": "PIP_G3", "G4": "PIP_G4"})
    gen_plot = gen_loci.merge(gen_wide, on="meta_locus_id", how="left")
    for c in ["PIP_G1", "PIP_G2", "PIP_G3", "PIP_G4"]:
        gen_plot[c] = pd.to_numeric(gen_plot.get(c, np.nan), errors="coerce").fillna(0)
    gen_plot["max_any_PIP"] = gen_plot[["PIP_G1", "PIP_G2", "PIP_G3", "PIP_G4"]].max(axis=1)
    gen_plot["locus_label"] = make_locus_label(gen_plot)
    gen_plot = gen_plot.sort_values(["max_any_PIP", "meta_locus_id"]).reset_index(drop=True)

    plt.figure(figsize=(8.5, max(4.5, 0.55 * len(gen_plot))))
    y = np.arange(len(gen_plot))
    for i, r in gen_plot.iterrows():
        vals = [r["PIP_G1"], r["PIP_G2"], r["PIP_G3"], r["PIP_G4"]]
        plt.plot([min(vals), max(vals)], [i, i], linewidth=1, alpha=0.6)
    plt.scatter(gen_plot["PIP_G1"], y, label="G1 Silent", s=34)
    plt.scatter(gen_plot["PIP_G2"], y, label="G2 Baby Boomer", s=34)
    plt.scatter(gen_plot["PIP_G3"], y, label="G3 Gen X", s=34)
    plt.scatter(gen_plot["PIP_G4"], y, label="G4 Millennial", s=34)
    plt.axvline(args.pip_threshold, linestyle=":", linewidth=1)
    plt.yticks(y, gen_plot["locus_label"], fontsize=8)
    plt.xlabel("Maximum PIP within core meta-locus")
    plt.ylabel("Generation-moderated meta-loci")
    plt.title("Generation-specific max-PIP patterns across\ngeneration-moderated meta-loci")
    plt.xlim(-0.02, 1.03)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "generation_moderated_maxPIP_core_plot.png"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "generation_moderated_maxPIP_core_plot.pdf"))
    plt.close()
    gen_plot.to_csv(os.path.join(args.out_dir, "generation_moderated_maxPIP_core_plot_data.tsv"), sep="\t", index=False)

    print("Saved plots to:", args.out_dir)


if __name__ == "__main__":
    main()
