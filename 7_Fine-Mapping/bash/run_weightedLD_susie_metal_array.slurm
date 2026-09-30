#!/usr/bin/env bash

set -euo pipefail


PROJECT=${path}
PGEN_DIR=${path}
ANALYSIS_ID=${path}
METAL_MODE=${path}
JOB_GROUP=${JOB_GROUP:}
RSCRIPT=${RSCRIPT:}
WEIGHTS_FILE=${WEIGHTS_FILE:}
OFFSET=${OFFSET:-0}
L=${L:-10}
ANCESTRIES=${ANCESTRIES:-"AFR AMR EAS EUR"}

OUT_BASE=${PROJECT}/outputs/${ANALYSIS_ID}/${METAL_MODE}
JOBLIST_DIR=${OUT_BASE}/joblists

if [ "${JOB_GROUP}" == "A" ]; then
  JOBLIST=${JOBLIST_DIR}/jobs_A_nsnps_le5000_minimal.tsv
  RESULT_DIR=${OUT_BASE}/susie_results_A
  TMP_LD_DIR=${OUT_BASE}/ld_temp_A
  MAC_FILTER=${MAC_FILTER:-2}
elif [ "${JOB_GROUP}" == "B" ]; then
  JOBLIST=${JOBLIST_DIR}/jobs_B_nsnps_5001_10000_minimal.tsv
  RESULT_DIR=${OUT_BASE}/susie_results_B
  TMP_LD_DIR=${OUT_BASE}/ld_temp_B
  MAC_FILTER=${MAC_FILTER:-5}
elif [ "${JOB_GROUP}" == "C_split" ]; then
  JOBLIST=${JOBLIST_DIR}/jobs_C_split_minimal.tsv
  RESULT_DIR=${OUT_BASE}/susie_results_C_split
  TMP_LD_DIR=${OUT_BASE}/ld_temp_C_split
  MAC_FILTER=${MAC_FILTER:-5}
else
  echo "Unknown JOB_GROUP=${JOB_GROUP}; expected A, B, or C_split"
  exit 1
fi

JOBLIST=${JOBLIST_OVERRIDE:-$JOBLIST}
echo "Using JOBLIST=${JOBLIST}"

mkdir -p logs "${RESULT_DIR}" "${TMP_LD_DIR}"

JOB_INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))
LINE=$(awk -v task_id="${JOB_INDEX}" 'BEGIN{FS=OFS="\t"} NR==task_id+1 {print}' "${JOBLIST}")

if [ -z "${LINE}" ]; then
  echo "No job found for JOB_INDEX=${JOB_INDEX}"
  exit 1
fi

if [ "${JOB_GROUP}" == "C_split" ]; then
  IFS=$'\t' read -r META_ANC ORIGINAL_LOCUS_ID LOCUS_ID CHR START END N_SNPS SUMSTATS_FILE SNPLIST_FILE <<< "${LINE}"
else
  IFS=$'\t' read -r META_ANC LOCUS_ID CHR N_SNPS SUMSTATS_FILE SNPLIST_FILE <<< "${LINE}"
  ORIGINAL_LOCUS_ID=${LOCUS_ID}
fi

echo "ANALYSIS_ID=${ANALYSIS_ID}"
echo "JOB_GROUP=${JOB_GROUP}"
echo "JOB_INDEX=${JOB_INDEX}"
echo "LOCUS_ID=${LOCUS_ID}"
echo "CHR=${CHR}"
echo "N_SNPS=${N_SNPS}"
echo "MAC_FILTER=${MAC_FILTER}"
echo "SUMSTATS_FILE=${SUMSTATS_FILE}"
echo "SNPLIST_FILE=${SNPLIST_FILE}"
echo "HOSTNAME=$(hostname)"

if ! [[ "${CHR}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: CHR is not numeric. LINE=${LINE}"
  exit 1
fi

META_TMP_DIR=${TMP_LD_DIR}/META
META_RESULT_DIR=${RESULT_DIR}/META
mkdir -p "${META_TMP_DIR}" "${META_RESULT_DIR}"

OUT_PREFIX=${META_RESULT_DIR}/${LOCUS_ID}_META

if [ -s "${OUT_PREFIX}.pip.tsv" ] && [ -s "${OUT_PREFIX}.summary.tsv" ]; then
  echo "Weighted-LD SuSiE result already exists. Skipping."
  exit 0
fi

LD_SPECS=()

for ANC in ${ANCESTRIES}; do
  REF=${PGEN_DIR}/1000G_GRCh38_${ANC}_chr${CHR}
  LD_PREFIX=${META_TMP_DIR}/${LOCUS_ID}_${ANC}

  echo "Generating LD for ${ANC}: ${REF}"

  plink2 \
    --pfile "${REF}" \
    --extract "${SNPLIST_FILE}" \
    --mac "${MAC_FILTER}" \
    --geno 0.05 \
    --r-unphased square zs ref-based \
    --out "${LD_PREFIX}"

  LD_ZST=${LD_PREFIX}.unphased.vcor1.zst
  LD_VARS=${LD_PREFIX}.unphased.vcor1.vars

  if [ ! -s "${LD_ZST}" ] || [ ! -s "${LD_VARS}" ]; then
    echo "WARNING: Missing LD outputs for ${ANC}; skipping this ancestry in weighted LD."
    continue
  fi

  LD_SPECS+=("${ANC}:${LD_ZST}:${LD_VARS}")
done

if [ "${#LD_SPECS[@]}" -eq 0 ]; then
  echo "ERROR: No ancestry LD matrices were generated."
  exit 1
fi

Rscript "${RSCRIPT}" \
  "${SUMSTATS_FILE}" \
  "${OUT_PREFIX}" \
  "${L}" \
  "${WEIGHTS_FILE}" \
  "${ANALYSIS_ID}" \
  "${LD_SPECS[@]}"

rm -f ${META_TMP_DIR}/${LOCUS_ID}_*.unphased.vcor1.zst
rm -f ${META_TMP_DIR}/${LOCUS_ID}_*.unphased.vcor1.vars

echo "Done: ${ANALYSIS_ID} ${JOB_GROUP} ${LOCUS_ID}"
