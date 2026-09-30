#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --job-name=popcorn_step1_hg19
#SBATCH --mem=40000
#SBATCH --output=/YOUR PATH HERE/Popcorn/popcorn_step1_hg19.out
#SBATCH --error=/YOUR PATH HERE/Popcorn/popcorn_step1_hg19.error
#SBATCH --mail-user=jeong212@purdue.edu
#SBATCH --mail-type=END,FAIL

module load anaconda/2022.10-py39
#source ~/envs/popcorn/bin/activate

#popcorn compute -h
#popcorn fit -h

soft_dir="/YOUR PATH HERE/Popcorn"
ref="hg19"
ref_dir="/YOUR PATH HERE/1000G/${ref}/clean"
root="/YOUR PATH HERE/"

#python -m popcorn compute -v 2 --bfile1 ${dir}/test/EUR_ref --bfile2 ${dir}/test/EAS_ref --gen_effect ${root}/Output/EUR_EAS_test_ge.cscore

declare -A anc_name

anc_name[EUR]="european"
anc_name[EAS]="east_asian"
anc_name[AFR]="african_american"
anc_name[AMR]="latino"

ancestries=(EUR EAS AMR AFR)


for ((i=0; i<${#ancestries[@]}; i++)); do
  for ((j=i+1; j<${#ancestries[@]}; j++)); do

    ANC1=${ancestries[$i]}
    ANC2=${ancestries[$j]}

    NAME1=${anc_name[$ANC1]}
    NAME2=${anc_name[$ANC2]}
    
    echo ""
    echo "*=== Running ${ANC1}_${ANC2} ===*"
    echo ""
    python ${soft_dir}/popcorn/__main__.py compute -v 2 \
	    --SNPs_to_store 500000 \
	    --bfile1 ${ref_dir}/${ANC1}/${ANC1}.chr1_22 \
	    --bfile2 ${ref_dir}/${ANC2}/${ANC2}.chr1_22 \
	    --gen_effect ${root}/Output/Step1_Scores/${ref}/${NAME1}_${NAME2}.cscore

    echo ""
    echo "*===========================*"
    echo ""

  done #for: j
done #for: i

