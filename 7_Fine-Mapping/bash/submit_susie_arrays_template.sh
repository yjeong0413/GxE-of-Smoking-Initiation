#!/usr/bin/env bash
set -euo pipefail

# Template for submitting SuSiE array jobs in chunks.
# Requires an external SLURM script and R script used for weighted-LD SuSiE.

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
META_BASE="${META_BASE:-${PROJECT_DIR}/outputs/meta_locus}"
PGEN_DIR="${PGEN_DIR:-/path/to/1000G_pgen_reference}"
METAL_MODE="${METAL_MODE:-meta_locus_400kb}"
SCRIPT="${SCRIPT:-/path/to/run_weightedLD_susie_metal_array.slurm}"
RSCRIPT="${RSCRIPT:-/path/to/run_susie_metal_weighted_ld_one_locus.R}"

ANALYSIS_ID="${ANALYSIS_ID:?Set ANALYSIS_ID}"
JOB_GROUP="${JOB_GROUP:?Set JOB_GROUP: A, B, or C_split}"
WEIGHTS_FILE="${WEIGHTS_FILE:?Set WEIGHTS_FILE}"
JOBLIST="${JOBLIST:?Set JOBLIST}"
CHUNK_SIZE="${CHUNK_SIZE:-500}"
MAX_CONCURRENT="${MAX_CONCURRENT:-50}"

N_JOBS=$(( $(wc -l < "${JOBLIST}") - 1 ))

for START in $(seq 1 "${CHUNK_SIZE}" "${N_JOBS}"); do
  END=$((START + CHUNK_SIZE - 1))
  if [ "${END}" -gt "${N_JOBS}" ]; then
    END=${N_JOBS}
  fi

  ARRAY_N=$((END - START + 1))
  OFFSET=$((START - 1))

  echo "Submitting ${ANALYSIS_ID} ${JOB_GROUP}: jobs ${START}-${END}"
  sbatch \
    --array=1-${ARRAY_N}%${MAX_CONCURRENT} \
    --export=ALL,PROJECT=${PROJECT_DIR},PGEN_DIR=${PGEN_DIR},ANALYSIS_ID=${ANALYSIS_ID},METAL_MODE=${METAL_MODE},JOB_GROUP=${JOB_GROUP},RSCRIPT=${RSCRIPT},WEIGHTS_FILE=${WEIGHTS_FILE},JOBLIST_OVERRIDE=${JOBLIST},OFFSET=${OFFSET} \
    "${SCRIPT}"
done
