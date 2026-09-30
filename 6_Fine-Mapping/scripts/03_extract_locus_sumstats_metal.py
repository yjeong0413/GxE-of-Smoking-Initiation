#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-id", required=True)
    parser.add_argument("--harm-dir", required=True)
    parser.add_argument("--union-file", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-snps", type=int, default=10)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    loci = pd.read_csv(args.union_file, sep="\t")
    loci.columns = loci.columns.str.strip()
    loci["CHR"] = loci["CHR"].astype(str).str.replace("^chr", "", regex=True)
    loci["START"] = pd.to_numeric(loci["START"], errors="raise").astype(int)
    loci["END"] = pd.to_numeric(loci["END"], errors="raise").astype(int)

    out_meta = os.path.join(args.out_dir, "META")
    snp_dir = os.path.join(out_meta, "snplists")
    os.makedirs(out_meta, exist_ok=True)
    os.makedirs(snp_dir, exist_ok=True)

    rows = []

    for chrom in range(1, 23):
        chrom = str(chrom)
        gwas_file = os.path.join(args.harm_dir, f"smoking_META_chr{chrom}_harmonized.tsv")

        chr_loci = loci[loci["CHR"].astype(str) == chrom].copy()

        if not os.path.exists(gwas_file):
            for _, locus in chr_loci.iterrows():
                rows.append({
                    "analysis_id": args.analysis_id,
                    "ancestry": "META",
                    "locus_id": locus["locus_id"],
                    "CHR": chrom,
                    "START": int(locus["START"]),
                    "END": int(locus["END"]),
                    "locus_size_bp": int(locus.get("locus_size_bp", int(locus["END"]) - int(locus["START"]) + 1)),
                    "n_snps": 0,
                    "status": "missing_chr_file",
                    "sumstats_file": "",
                    "snplist_file": ""
                })
            continue

        gwas = pd.read_csv(gwas_file, sep="\t")
        gwas["BP"] = pd.to_numeric(gwas["BP"], errors="coerce")
        gwas = gwas.dropna(subset=["BP"]).copy()
        gwas["BP"] = gwas["BP"].astype(int)

        for _, locus in chr_loci.iterrows():
            locus_id = locus["locus_id"]
            start = int(locus["START"])
            end = int(locus["END"])

            d = gwas[(gwas["BP"] >= start) & (gwas["BP"] <= end)].copy()
            n_snps = len(d)

            status = "kept" if n_snps >= args.min_snps else "too_few_snps"

            sumstats_file = os.path.join(out_meta, f"{locus_id}_META.sumstats.tsv")
            snplist_file = os.path.join(snp_dir, f"{locus_id}_META.snplist.txt")

            if status == "kept":
                d = d.sort_values("BP").drop_duplicates("variant_id", keep="first")
                d.to_csv(sumstats_file, sep="\t", index=False)
                d["variant_id"].to_csv(snplist_file, index=False, header=False)

            rows.append({
                "analysis_id": args.analysis_id,
                "ancestry": "META",
                "locus_id": locus_id,
                "CHR": chrom,
                "START": start,
                "END": end,
                "locus_size_bp": int(locus.get("locus_size_bp", end - start + 1)),
                "n_snps": n_snps,
                "status": status,
                "sumstats_file": sumstats_file if status == "kept" else "",
                "snplist_file": snplist_file if status == "kept" else ""
            })

    summary = pd.DataFrame(rows)
    out_summary = os.path.join(args.out_dir, "locus_sumstats_extraction_summary.tsv")
    summary.to_csv(out_summary, sep="\t", index=False)

    print("Saved:", out_summary)
    print(summary["status"].value_counts())
    print(summary["n_snps"].describe())

if __name__ == "__main__":
    main()
