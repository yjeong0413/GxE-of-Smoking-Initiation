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

echo " "
# 2) Define directories
echo "Define directories."
root_dir=/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc/${ANCESTRY}
dest_dir=${root_dir}/CLEAN/GWAS_QC
output_raw_dir=${dest_dir}/Genetic_Correlation/raw/${PHENOTYPE}
output_liability_dir=${dest_dir}//Genetic_Correlation/liability/${PHENOTYPE}

echo "Done."

zcat  $dest_dir/ever_tobacco_user_female/ever_tobacco_user_female.sumstats.gz | head -5

# Raw - rg(female, male)
subgroup="gender"
SUBGROUP1="female"
SUBGROUP2="male"

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_${SUBGROUP1}/ever_tobacco_user_${SUBGROUP1}.sumstats.gz,$dest_dir/ever_tobacco_user_${SUBGROUP2}/ever_tobacco_user_${SUBGROUP2}.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_${SUBGROUP1}_${SUBGROUP2}_rg_raw


subgroup="region"

# Raw - rg(MW,NE), rg(MW,SE), rg(MW,SW), rg(MW,W)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_MW/ever_tobacco_user_MW.sumstats.gz,$dest_dir/ever_tobacco_user_NE/ever_tobacco_user_NE.sumstats.gz,$dest_dir/ever_tobacco_user_SE/ever_tobacco_user_SE.sumstats.gz,$dest_dir/ever_tobacco_user_SW/ever_tobacco_user_SW.sumstats.gz,$dest_dir/ever_tobacco_user_W/ever_tobacco_user_W.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_MWNESESWW_rg_raw

# Raw - rg(NE,SE), rg(NE,SW), rg(NE,W)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_NE/ever_tobacco_user_NE.sumstats.gz,$dest_dir/ever_tobacco_user_SE/ever_tobacco_user_SE.sumstats.gz,$dest_dir/ever_tobacco_user_SW/ever_tobacco_user_SW.sumstats.gz,$dest_dir/ever_tobacco_user_W/ever_tobacco_user_W.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_NESESWW_rg_raw

# Raw - rg(SE,SW), rg(SE,W)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_SE/ever_tobacco_user_SE.sumstats.gz,$dest_dir/ever_tobacco_user_SW/ever_tobacco_user_SW.sumstats.gz,$dest_dir/ever_tobacco_user_W/ever_tobacco_user_W.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_SESWW_rg_raw

# Raw - rg(SW,W)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_SW/ever_tobacco_user_SW.sumstats.gz,$dest_dir/ever_tobacco_user_W/ever_tobacco_user_W.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_SWW_rg_raw

subgroup="birthyear"

# Raw - rg(bin4,bin5),..., rg(bin4,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin4/ever_tobacco_user_birthyear_bin4.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin5/ever_tobacco_user_birthyear_bin5.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin6/ever_tobacco_user_birthyear_bin6.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin7/ever_tobacco_user_birthyear_bin7.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin8/ever_tobacco_user_birthyear_bin8.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin419_rg_raw

# Raw - rg(bin5,bin6),..., rg(bin5,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin5/ever_tobacco_user_birthyear_bin5.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin6/ever_tobacco_user_birthyear_bin6.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin7/ever_tobacco_user_birthyear_bin7.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin8/ever_tobacco_user_birthyear_bin8.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin519_rg_raw

# Raw - rg(bin6,bin7),..., rg(bin6,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin6/ever_tobacco_user_birthyear_bin6.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin7/ever_tobacco_user_birthyear_bin7.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin8/ever_tobacco_user_birthyear_bin8.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin619_rg_raw

# Raw - rg(bin7,bin8),..., rg(bin7,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin7/ever_tobacco_user_birthyear_bin7.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin8/ever_tobacco_user_birthyear_bin8.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin719_rg_raw

# Raw - rg(bin8,bin9),..., rg(bin8,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin8/ever_tobacco_user_birthyear_bin8.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin819_rg_raw

# Raw - rg(bin9,bin10),..., rg(bin9,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin9/ever_tobacco_user_birthyear_bin9.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin919_rg_raw

# Raw - rg(bin10,bin11),..., rg(bin10,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin10/ever_tobacco_user_birthyear_bin10.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1019_rg_raw

# Raw - rg(bin11,bin12),..., rg(bin11,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin11/ever_tobacco_user_birthyear_bin11.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1119_rg_raw

# Raw - rg(bin12,bin13),..., rg(bin12,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin12/ever_tobacco_user_birthyear_bin12.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1219_rg_raw

# Raw - rg(bin13,bin14),..., rg(bin13,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin13/ever_tobacco_user_birthyear_bin13.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1319_rg_raw

# Raw - rg(bin14,bin15),..., rg(bin14,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin14/ever_tobacco_user_birthyear_bin14.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1419_rg_raw

# Raw - rg(bin15,bin16),..., rg(bin15,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin15/ever_tobacco_user_birthyear_bin15.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1519_rg_raw

# Raw - rg(bin16,bin17),..., rg(bin16,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin16/ever_tobacco_user_birthyear_bin16.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1619_rg_raw

# Raw - rg(bin17,bin18),..., rg(bin17,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin17/ever_tobacco_user_birthyear_bin17.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1719_rg_raw

# Raw - rg(bin18,bin19)

ldsc.py \
      --rg $dest_dir/ever_tobacco_user_birthyear_bin18/ever_tobacco_user_birthyear_bin18.sumstats.gz,$dest_dir/ever_tobacco_user_birthyear_bin19/ever_tobacco_user_birthyear_bin19.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_bin1819_rg_raw

zcat  $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN4.sumstats.gz | head -5

# Raw - rg(bin1,bin2), rg(bin1,bin3), rg(bin1,bin4)

ldsc.py \
      --rg $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN1.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN2.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN3.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN4.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_groupingBIN14_rg_raw

# Raw - rg(bin2,bin3), rg(bin2,bin4)

ldsc.py \
      --rg $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN2.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN3.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN4.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_groupingBIN24_rg_raw

# Raw - rg(bin3,bin4)

ldsc.py \
      --rg $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN3.sumstats.gz,$dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN4.sumstats.gz \
      --ref-ld-chr $root_dir/eur_w_ld_chr/ \
      --w-ld-chr $root_dir/eur_w_ld_chr/ \
      --out $output_raw_dir/ever_tobacco_user_${subgroup}/ever_tobacco_user_birthyear_groupingBIN34_rg_raw


