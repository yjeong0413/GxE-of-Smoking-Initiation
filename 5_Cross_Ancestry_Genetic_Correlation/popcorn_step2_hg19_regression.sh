#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --job-name=popcorn_step2_hg19_regression
#SBATCH --mem=40000
#SBATCH --output=/YOUR PATH HERE/Popcorn/popcorn_step2_hg19_regression.out
#SBATCH --error=/YOUR PATH HERE/Popcorn/popcorn_step2_hg19_regression.error
#SBATCH --mail-user=jeong212@purdue.edu
#SBATCH --mail-type=END,FAIL

module load anaconda/2022.10-py39
#source ~/envs/popcorn/bin/activate

#popcorn compute -h
#popcorn fit -h

soft_dir="/YOUR PATH HERE/Popcorn"
root="/YOUR PATH HERE/Popcorn"
ref="hg19"

declare -A anc_name

anc_name[EUR]="european"
anc_name[EAS]="east_asian"
anc_name[AFR]="african_american"
anc_name[AMR]="latino"

ancestries=(EUR EAS AMR AFR)

#[GENDER]#
echo ""
echo "#######################[GENDER]########################"
echo ""
type="gender"
subgroups=(female male)

for SUBGROUP in "${subgroups[@]}"; do
  echo ""
  echo "*=== ${SUBGROUP} ===*"
  echo ""
for ((i=0; i<${#ancestries[@]}; i++)); do
  for ((j=i+1; j<${#ancestries[@]}; j++)); do

    ANC1=${ancestries[$i]}
    ANC2=${ancestries[$j]}

    NAME1=${anc_name[$ANC1]}
    NAME2=${anc_name[$ANC2]}

    echo ""
    echo "*--- Running ${ANC1}_${ANC2} ---*"
    echo ""
    python ${soft_dir}/popcorn/__main__.py fit -v 1 \
		--cfile ${root}/Output//Step1_Scores/${ref}/${NAME1}_${NAME2}_cleaned.cscore \
		--gen_effect \
		--sfile1 ${root}/Sumstats/${NAME1}/${type}/${NAME1}_${SUBGROUP}.txt \
		--sfile2 ${root}/Sumstats/${NAME2}/${type}/${NAME2}_${SUBGROUP}.txt \
		${root}/Output/Step2_Results/${ref}/${SUBGROUP}/${NAME1}_${NAME2}_ge_results_${SUBGROUP}_regression.txt
    echo ""
    echo "*--------------------------------*"
    echo ""

  done #for: j
done #for: i
  echo ""
  echo "=============================="
  echo ""
done #for: SUBGROUP
echo ""
echo "#####################################################"
echo ""



#[GENERATION]#
echo ""
echo "####################[GENERATION]#####################"
echo ""
type="birthyear"
subgroups=(groupingBIN{1..4})

for SUBGROUP in "${subgroups[@]}"; do
  echo ""
  echo "*=== ${SUBGROUP} ===*"
  echo ""
for ((i=0; i<${#ancestries[@]}; i++)); do
  for ((j=i+1; j<${#ancestries[@]}; j++)); do

    ANC1=${ancestries[$i]}
    ANC2=${ancestries[$j]}

    NAME1=${anc_name[$ANC1]}
    NAME2=${anc_name[$ANC2]}

    echo ""
    echo "*--- Running ${ANC1}_${ANC2} ---*"
    echo ""
    python ${soft_dir}/popcorn/__main__.py fit -v 1 \
		--cfile ${root}/Output/Step1_Scores/${ref}/${NAME1}_${NAME2}_cleaned.cscore \
		--gen_effect \
		--sfile1 ${root}/Sumstats/${NAME1}/${type}/${NAME1}_${SUBGROUP}.txt \
		--sfile2 ${root}/Sumstats/${NAME2}/${type}/${NAME2}_${SUBGROUP}.txt \
		${root}/Output/Step2_Results/${ref}/${SUBGROUP}/${NAME1}_${NAME2}_ge_results_${SUBGROUP}_regression.txt
    echo ""
    echo "*--------------------------------*"
    echo ""

  done #for: j
done #for: i
  echo ""
  echo "=============================="
  echo ""
done #for: SUBGROUP
echo ""
echo "######################################################"
echo ""
