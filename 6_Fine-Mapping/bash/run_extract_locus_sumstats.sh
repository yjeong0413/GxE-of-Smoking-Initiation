#!/usr/bin/env bash
set -euo pipefail

# Template: extract locus-level summary statistics for each analysis.
# Edit paths or pass environment variables from the command line.

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
META_BASE="${META_BASE:-${PROJECT_DIR}/outputs/meta_locus}"
UNION_FILE="${UNION_FILE:-${META_BASE}/meta_loci_union_hg38_buffered_for_sumstats.tsv}"
EXTRACT_SCRIPT="${EXTRACT_SCRIPT:-${PROJECT_DIR}/scripts/03_extract_locus_sumstats_metal.py}"

ANALYSES=(
  sex_female_meta
  sex_male_meta
  generation1_meta
  generation2_meta
  generation3_meta
  generation4_meta
)

for ANALYSIS_ID in "${ANALYSES[@]}"; do
  echo "Extracting ${ANALYSIS_ID}"
  HARM_DIR="${PROJECT_DIR}/outputs/${ANALYSIS_ID}/harmonized_META_N80"
  OUT_DIR="${META_BASE}/${ANALYSIS_ID}/locus_sumstats"

  python "${EXTRACT_SCRIPT}" \
    --analysis-id "${ANALYSIS_ID}" \
    --harm-dir "${HARM_DIR}" \
    --union-file "${UNION_FILE}" \
    --out-dir "${OUT_DIR}" \
    --min-snps 10
done
