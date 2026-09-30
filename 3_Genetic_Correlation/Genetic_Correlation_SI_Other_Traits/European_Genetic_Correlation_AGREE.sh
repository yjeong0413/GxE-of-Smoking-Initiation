#Load the necessary modules (adjust based on your HPC system)
module load biocontainers/default
echo "Loading [biocontainers/default]: Done."
module load ldsc
echo "loading [ldsc]: done."

echo " "
echo "Check whether [ldsc] is well loaded."
ldsc.py -h
echo " "
echo "Check whether [munge_sumstats] is well loaded."
munge_sumstats.py -h
echo "Done."

# 1) Define Ancestry and Phenotype
ANCESTRY="european"
PHENOTYPE="SI"
OTHER="AGREE"
GWAS_LIST=("Gupta_et_al_2024")

echo " "
# 2) Define directories
echo "Define directories."
root_dir=/YOUR PATH HERE/
dest_dir=${root_dir}/${ANCESTRY}/CLEAN/GWAS_QC
munge_stats_Other_dir=${root_dir}/Other_Sumstats
output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
#output_liability_dir=${dest_dir}//Genetic_Correlation/liability/${PHENOTYPE}

echo "Done."

# Ancestry Setting
ref_ancestry=$root_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr  # European

# Phenotype Setting
munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}

# Check first three lines of [ever_tobacco_user.sumstats.gz]
zcat ${dest_dir}/ever_tobacco_user/ever_tobacco_user.sumstats.gz | head -n 4

# Check first three lines of [${OTHER}_${GWAS}]
for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}
    
    echo "[${OTHER}_${GWAS}.sumstats.gz]"
    zcat $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz | head -n 4
    echo
done  #for: GWAS

subgroup="full"

for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}

 ldsc.py \
        --rg $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz,${dest_dir}/ever_tobacco_user/ever_tobacco_user.sumstats.gz \
        --ref-ld-chr $ref_ancestry/ \
        --w-ld-chr $ref_ancestry/ \
        --out $output_raw_dir/${subgroup}/ever_tobacco_user_${subgroup}_${OTHER}_rg_raw
done #for: GWAS

# Check first three lines of [ever_tobacco_user_female.sumstats.gz]
echo "Female"
zcat ${dest_dir}/ever_tobacco_user_female/ever_tobacco_user_female.sumstats.gz | head -n 4

# Check first three lines of [ever_tobacco_user_male.sumstats.gz]
echo "Male"
zcat ${dest_dir}/ever_tobacco_user_male/ever_tobacco_user_male.sumstats.gz | head -n 4

for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}
    
    echo "[${OTHER}_${GWAS}.sumstats.gz]"
    zcat $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz | head -n 4
    echo
done  #for: GWAS

subgroup="gender"
SUBGROUP1="female"
SUBGROUP2="male"

for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}

 ldsc.py \
        --rg $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz,${dest_dir}/ever_tobacco_user_${SUBGROUP1}/ever_tobacco_user_${SUBGROUP1}.sumstats.gz,${dest_dir}/ever_tobacco_user_${SUBGROUP2}/ever_tobacco_user_${SUBGROUP2}.sumstats.gz \
        --ref-ld-chr $ref_ancestry/ \
        --w-ld-chr $ref_ancestry/ \
        --out $output_raw_dir/${subgroup}/ever_tobacco_user_${subgroup}_${OTHER}_rg_raw
done #for: GWAS

# Check first three lines of [ever_tobacco_user_MW.sumstats.gz]
echo "MW"
zcat ${dest_dir}/ever_tobacco_user_MW/ever_tobacco_user_MW.sumstats.gz | head -n 4

# Check first three lines of [ever_tobacco_user_NE.sumstats.gz]
echo "NE"
zcat ${dest_dir}/ever_tobacco_user_NE/ever_tobacco_user_NE.sumstats.gz | head -n 4

# Check first three lines of [ever_tobacco_user_SE.sumstats.gz]
echo "SE"
zcat ${dest_dir}/ever_tobacco_user_SE/ever_tobacco_user_SE.sumstats.gz | head -n 4

# Check first three lines of [ever_tobacco_user_SW.sumstats.gz]
echo "SW"
zcat ${dest_dir}/ever_tobacco_user_SW/ever_tobacco_user_SW.sumstats.gz | head -n 4

# Check first three lines of [ever_tobacco_user_W.sumstats.gz]
echo "W"
zcat ${dest_dir}/ever_tobacco_user_W/ever_tobacco_user_W.sumstats.gz | head -n 4

# Check first three lines of [EA_Okbay_et_al_2022.sumstats.gz]
for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}
    
    echo "[${OTHER}_${GWAS}.sumstats.gz]"
    zcat $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz | head -n 4
    echo
done  #for: GWAS

subgroup="region"
for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}

 ldsc.py \
        --rg ${munge_Other_dir}/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz,${dest_dir}/ever_tobacco_user_MW/ever_tobacco_user_MW.sumstats.gz,${dest_dir}/ever_tobacco_user_NE/ever_tobacco_user_NE.sumstats.gz,${dest_dir}/ever_tobacco_user_SE/ever_tobacco_user_SE.sumstats.gz,${dest_dir}/ever_tobacco_user_SW/ever_tobacco_user_SW.sumstats.gz,${dest_dir}/ever_tobacco_user_W/ever_tobacco_user_W.sumstats.gz \
        --ref-ld-chr $ref_ancestry/ \
        --w-ld-chr $ref_ancestry/ \
        --out $output_raw_dir/${subgroup}/ever_tobacco_user_${subgroup}_${OTHER}_rg_raw
done #for:GWAS

for SUBGROUP in {4..19}; do

    # Check first three lines of [ever_tobacco_user_{subgroup}.sumstats.gz]
    echo "bin${SUBGROUP}"
    zcat ${dest_dir}/ever_tobacco_user_birthyear_bin${SUBGROUP}/ever_tobacco_user_birthyear_bin${SUBGROUP}.sumstats.gz | head -n 4
done

# Check first three lines of [EA_Okbay_et_al_2022.sumstats.gz]
for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}
    
    echo "[${OTHER}_${GWAS}.sumstats.gz]"
    zcat $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz | head -n 4
    echo
done  #for: GWAS

subgroup="birthyear"

for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}

 ldsc.py \
        --rg $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin4/ever_tobacco_user_birthyear_bin4.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin5/ever_tobacco_user_birthyear_bin5.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin6/ever_tobacco_user_birthyear_bin6.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin7/ever_tobacco_user_birthyear_bin7.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin8/ever_tobacco_user_birthyear_bin8.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,${dest_dir}/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz  \
        --ref-ld-chr $ref_ancestry/ \
        --w-ld-chr $ref_ancestry/ \
        --out $output_raw_dir/${subgroup}/ever_tobacco_user_${subgroup}_${OTHER}_rg_raw
done #for: GWAS

for SUBGROUP in {1..4}; do
    # Check first three lines of [METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN${SUBGROUP}.sumstats.gz]
    echo "[groupingBin${SUBGROUP}]"
    zcat  $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN${SUBGROUP}.sumstats.gz | head -5
    echo ""
done

subgroup="birthyear"

for GWAS in "${GWAS_LIST[@]}"; do
    output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${OTHER}/${GWAS}
    munge_Other_dir=${munge_stats_Other_dir}/${OTHER}/${GWAS}

ldsc.py \
      --rg  $munge_Other_dir/munge_sumstats/${OTHER}_${GWAS}.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN1.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN2.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN3.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN4.sumstats.gz \
        --ref-ld-chr $ref_ancestry/ \
        --w-ld-chr $ref_ancestry/ \
      --out $output_raw_dir/${subgroup}/ever_tobacco_user_${subgroup}_${OTHER}_rg_raw_2ndMA
done #for: GWAS


