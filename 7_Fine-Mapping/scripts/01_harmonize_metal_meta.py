#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np

VALID = {"A", "C", "G", "T"}
AMBIG = {frozenset(["A", "T"]), frozenset(["C", "G"])}

def clean_chr(x):
    return (
        x.astype(str)
        .str.replace("^chr", "", regex=True)
        .str.replace("^CHR", "", regex=True)
    )

def choose_col(df, candidates, required=True):
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    if required:
        raise ValueError(f"None of these columns found: {candidates}. Available: {list(df.columns)}")
    return None

def parse_marker(marker_series):
    parts = marker_series.astype(str).str.split(":", expand=True)
    if parts.shape[1] < 2:
        raise ValueError("Marker column cannot be parsed as CHR:BP")
    return parts[0], parts[1]

def is_ambiguous(ref, alt):
    return frozenset([ref, alt]) in AMBIG

def find_map_file(map_dir, chrom):
    candidates = [
        os.path.join(map_dir, f"1000G_GRCh38_ALL_chr{chrom}.variant_map.tsv"),
        os.path.join(map_dir, f"1000G_GRCh38_EUR_chr{chrom}.variant_map.tsv"),
        os.path.join(map_dir, f"chr{chrom}.variant_map.tsv"),
    ]
    for f in candidates:
        if os.path.exists(f):
            return f
    raise FileNotFoundError(f"No variant map found for chr{chrom}. Tried: {candidates}")

def load_ref_map(map_dir, chrom):
    f = find_map_file(map_dir, chrom)
    ref = pd.read_csv(f, sep="\t")
    ref.columns = ref.columns.str.strip()

    required = ["CHR", "BP", "variant_id", "REF", "ALT"]
    missing = [c for c in required if c not in ref.columns]
    if missing:
        raise ValueError(f"Missing columns in map {f}: {missing}")

    ref = ref[required].copy()
    ref["CHR"] = clean_chr(ref["CHR"])
    ref["BP"] = pd.to_numeric(ref["BP"], errors="coerce").astype("Int64")
    ref["REF"] = ref["REF"].astype(str).str.upper()
    ref["ALT"] = ref["ALT"].astype(str).str.upper()
    ref = ref.dropna(subset=["BP"])
    return ref

def harmonize_one(analysis_id, gwas_file, map_dir, out_dir, chunksize):
    os.makedirs(out_dir, exist_ok=True)

    # remove previous outputs for this analysis
    for chrom in range(1, 23):
        f = os.path.join(out_dir, f"smoking_META_chr{chrom}_harmonized.tsv")
        if os.path.exists(f):
            os.remove(f)

    first = pd.read_csv(gwas_file, sep=r"\s+", engine="python", nrows=5)
    first.columns = first.columns.str.strip()

    chr_col = choose_col(first, ["CHR", "chr", "#CHROM", "Chromosome", "chromosome"], required=False)
    bp_col = choose_col(first, ["BP", "bp", "POS", "pos", "Position", "position"], required=False)
    marker_col = choose_col(first, ["MarkerName", "Marker", "SNP", "snpid", "variant_id", "ID"], required=False)

    ea_col = choose_col(first, ["Allele1", "A1", "EA", "effect_allele", "B"])
    oa_col = choose_col(first, ["Allele2", "A2", "OA", "other_allele", "A"])
    beta_col = choose_col(first, ["Effect_new", "Effect", "BETA", "Beta", "beta", "effect"])
    se_col = choose_col(first, ["StdErr", "SE", "stderr", "se"])
    p_col = choose_col(first, ["P-value", "P", "PVAL", "pval", "p_value"])
    n_col = choose_col(first, ["N", "N_tot", "TotalSampleSize", "N_total", "n"], required=True)

    if chr_col is None or bp_col is None:
        if marker_col is None:
            raise ValueError("Need CHR/BP columns or marker column with CHR:BP.")
        print(f"CHR/BP not found; parsing from marker column: {marker_col}")

    usecols = list(set([x for x in [chr_col, bp_col, marker_col, ea_col, oa_col, beta_col, se_col, p_col, n_col] if x is not None]))

    summary = {
        str(chrom): {
            "analysis_id": analysis_id,
            "CHR": str(chrom),
            "input_rows_chr": 0,
            "after_basic_qc": 0,
            "position_matched": 0,
            "allele_matched": 0,
            "after_ambiguous_removed": 0,
            "final_output_rows": 0,
            "output_file": os.path.join(out_dir, f"smoking_META_chr{chrom}_harmonized.tsv")
        }
        for chrom in range(1, 23)
    }

    header_written = {str(chrom): False for chrom in range(1, 23)}
    ref_cache = {}

    reader = pd.read_csv(gwas_file, sep=r"\s+", engine="python", usecols=usecols, chunksize=chunksize)

    for i, chunk in enumerate(reader, start=1):
        print(f"{analysis_id}: chunk {i}")

        chunk.columns = chunk.columns.str.strip()

        if chr_col is not None and bp_col is not None:
            chunk["CHR"] = clean_chr(chunk[chr_col])
            chunk["BP"] = pd.to_numeric(chunk[bp_col], errors="coerce").astype("Int64")
        else:
            c, b = parse_marker(chunk[marker_col])
            chunk["CHR"] = clean_chr(c)
            chunk["BP"] = pd.to_numeric(b, errors="coerce").astype("Int64")

        chunk["EA"] = chunk[ea_col].astype(str).str.upper()
        chunk["OA"] = chunk[oa_col].astype(str).str.upper()
        chunk["BETA"] = pd.to_numeric(chunk[beta_col], errors="coerce")
        chunk["SE"] = pd.to_numeric(chunk[se_col], errors="coerce")
        chunk["P"] = pd.to_numeric(chunk[p_col], errors="coerce")
        chunk["N"] = pd.to_numeric(chunk[n_col], errors="coerce")
        chunk["Z_EA"] = chunk["BETA"] / chunk["SE"]

        if marker_col is not None:
            chunk["snpid"] = chunk[marker_col].astype(str)
        else:
            chunk["snpid"] = chunk["CHR"].astype(str) + ":" + chunk["BP"].astype(str)

        for chrom, gwas_chr in chunk.groupby("CHR"):
            chrom = str(chrom).replace("chr", "")
            if chrom not in summary:
                continue

            summary[chrom]["input_rows_chr"] += len(gwas_chr)

            gwas_chr = gwas_chr[
                gwas_chr["EA"].isin(VALID)
                & gwas_chr["OA"].isin(VALID)
                & gwas_chr["BP"].notna()
                & gwas_chr["BETA"].notna()
                & gwas_chr["SE"].notna()
                & gwas_chr["P"].notna()
                & gwas_chr["N"].notna()
                & (gwas_chr["SE"] > 0)
            ].copy()

            summary[chrom]["after_basic_qc"] += len(gwas_chr)

            if gwas_chr.empty:
                continue

            if chrom not in ref_cache:
                ref_cache[chrom] = load_ref_map(map_dir, chrom)

            merged = gwas_chr.merge(ref_cache[chrom], on=["CHR", "BP"], how="inner")
            summary[chrom]["position_matched"] += len(merged)

            if merged.empty:
                continue

            same = (merged["EA"] == merged["REF"]) & (merged["OA"] == merged["ALT"])
            flip = (merged["EA"] == merged["ALT"]) & (merged["OA"] == merged["REF"])

            merged["allele_status"] = "mismatch"
            merged.loc[same, "allele_status"] = "same"
            merged.loc[flip, "allele_status"] = "flip"

            matched = merged[merged["allele_status"].isin(["same", "flip"])].copy()
            summary[chrom]["allele_matched"] += len(matched)

            if matched.empty:
                continue

            matched["BETA_REF"] = np.where(matched["allele_status"] == "same", matched["BETA"], -matched["BETA"])
            matched["Z_REF"] = np.where(matched["allele_status"] == "same", matched["Z_EA"], -matched["Z_EA"])

            matched["ambiguous"] = [is_ambiguous(r, a) for r, a in zip(matched["REF"], matched["ALT"])]
            matched = matched[~matched["ambiguous"]].copy()
            summary[chrom]["after_ambiguous_removed"] += len(matched)

            if matched.empty:
                continue

            matched = matched.sort_values("P").drop_duplicates("variant_id", keep="first")

            keep_cols = [
                "variant_id", "snpid", "CHR", "BP", "REF", "ALT",
                "EA", "OA", "BETA", "SE", "BETA_REF", "Z_REF",
                "P", "N", "allele_status"
            ]

            out_file = summary[chrom]["output_file"]

            matched[keep_cols].to_csv(
                out_file,
                sep="\t",
                index=False,
                mode="a",
                header=not header_written[chrom]
            )

            header_written[chrom] = True

    for chrom in range(1, 23):
        chrom = str(chrom)
        f = summary[chrom]["output_file"]
        if os.path.exists(f):
            summary[chrom]["final_output_rows"] = sum(1 for _ in open(f)) - 1

    summary_df = pd.DataFrame(summary.values())
    summary_file = os.path.join(out_dir, "harmonization_summary_META.tsv")
    summary_df.to_csv(summary_file, sep="\t", index=False)

    print("Saved:", summary_file)
    print("Total final rows:", summary_df["final_output_rows"].sum())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--map-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--chunksize", type=int, default=100000)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    manifest.columns = manifest.columns.str.strip()
    manifest["analysis_id"] = manifest["analysis_id"].astype(str).str.strip()
    manifest["gwas_file"] = manifest["gwas_file"].astype(str).str.strip()

    row = manifest[manifest["analysis_id"] == args.analysis_id]
    if row.empty:
        raise ValueError(f"No manifest row for {args.analysis_id}")

    gwas_file = row.iloc[0]["gwas_file"]
    if not os.path.exists(gwas_file):
        raise FileNotFoundError(gwas_file)

    harmonize_one(args.analysis_id, gwas_file, args.map_dir, args.out_dir, args.chunksize)

if __name__ == "__main__":
    main()
