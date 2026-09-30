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

# 1) Define Ancestry
ANCESTRY="european"

echo " "
# 2) Define directories
echo "Define directories."
root_dir=/YOUR PATH HERE/${ANCESTRY}
dest_dir=${root_dir}/CLEAN/GWAS_QC
output_raw_dir=${dest_dir}/Heritability/raw
output_liability_dir=${dest_dir}//Heritability/liability

echo "Done."

# Raw
ldsc.py \
        --h2 $dest_dir/ever_tobacco_user/ever_tobacco_user.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --w-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --out $output_raw_dir/ever_tobacco_user_h2_raw

# Raw
ldsc.py \
        --h2 $dest_dir/ever_tobacco_user/ever_tobacco_user.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/baseline/baseline. \
        --frqfile-chr $root_dir/Source_Files/1000G_frq/1000G.mac5eur. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_h2_partitioned_finucane2015_raw

# Raw
ldsc.py \
        --h2 $dest_dir/ever_tobacco_user/ever_tobacco_user.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --frqfile-chr $root_dir/Source_Files/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_h2_partitioned_gazal2018_raw

cd $root_dir # This is for [--ref-ld-chr-cts $root_dir/Multi_tissue_gene_expr.ldcts]. Otherwise, we face a path problem.

# Raw
ldsc.py \
        --h2-cts $dest_dir/ever_tobacco_user/ever_tobacco_user.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --ref-ld-chr-cts $root_dir/Source_Files/Multi_tissue_gene_expr.ldcts \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_multi_tissue_gene_expr_raw

#Generation
subsample=(female male)

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Gender: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --w-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Gender: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/baseline/baseline. \
        --frqfile-chr $root_dir/Source_Files/1000G_frq/1000G.mac5eur. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_partitioned_finucane2015_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Gender: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --frqfile-chr $root_dir/Source_Files/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_partitioned_gazal2018_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

cd $root_dir # This is for [--ref-ld-chr-cts $root_dir/Multi_tissue_gene_expr.ldcts]. Otherwise, we face a path problem.

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Gender: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2-cts $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --ref-ld-chr-cts $root_dir/Source_Files/Multi_tissue_gene_expr.ldcts \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_multi_tissue_gene_expr_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

#Generation
subsample=(MW NE SE SW W)

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Region: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --w-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Region: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/baseline/baseline. \
        --frqfile-chr $root_dir/Source_Files/1000G_frq/1000G.mac5eur. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_partitioned_finucane2015_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Region: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --frqfile-chr $root_dir/Source_Files/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_partitioned_gazal2018_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

cd $root_dir # This is for [--ref-ld-chr-cts $root_dir/Multi_tissue_gene_expr.ldcts]. Otherwise, we face a path problem.

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Region: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2-cts $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --ref-ld-chr-cts $root_dir/Source_Files/Multi_tissue_gene_expr.ldcts \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_multi_tissue_gene_expr_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

#Generation
subsample=(birthyear_bin{4..19})

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: bin${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --w-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: bin${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/baseline/baseline. \
        --frqfile-chr $root_dir/Source_Files/1000G_frq/1000G.mac5eur. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_partitioned_finucane2015_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: bin${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --frqfile-chr $root_dir/Source_Files/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_h2_partitioned_gazal2018_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

cd $root_dir # This is for [--ref-ld-chr-cts $root_dir/Multi_tissue_gene_expr.ldcts]. Otherwise, we face a path problem.

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: bin${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2-cts $dest_dir/ever_tobacco_user_${SUBGROUP}/ever_tobacco_user_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --ref-ld-chr-cts $root_dir/Source_Files/Multi_tissue_gene_expr.ldcts \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_${SUBGROUP}_multi_tissue_gene_expr_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

#Generation
subsample=(groupingBIN{1..4})

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --w-ld-chr $root_dir/Source_Files/eur_w_ld_chr/ \
        --out $output_raw_dir/ever_tobacco_user_birthyear_${SUBGROUP}_h2_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/baseline/baseline. \
        --frqfile-chr $root_dir/Source_Files/1000G_frq/1000G.mac5eur. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_birthyear_${SUBGROUP}_h2_partitioned_finucane2015_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: ${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2 $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${SUBGROUP}.sumstats.gz \
        --overlap-annot \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --frqfile-chr $root_dir/Source_Files/1000G_Phase3_frq/1000G.EUR.QC. \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_birthyear_${SUBGROUP}_h2_partitioned_gazal2018_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

cd $root_dir # This is for [--ref-ld-chr-cts $root_dir/Multi_tissue_gene_expr.ldcts]. Otherwise, we face a path problem.

# Raw
for SUBGROUP in "${subsample[@]}"
do
 echo "*------------------------------- Birthyear: bin${SUBGROUP} -------------------------------*"
 echo " "
 ldsc.py \
        --h2-cts $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${SUBGROUP}.sumstats.gz \
        --ref-ld-chr $root_dir/Source_Files/1000G_EUR_Phase3_baseline/baselineLD. \
        --ref-ld-chr-cts $root_dir/Source_Files/Multi_tissue_gene_expr.ldcts \
        --w-ld-chr $root_dir/Source_Files/weights_hm3_no_hla/weights. \
        --out $output_raw_dir/ever_tobacco_user_birthyear_${SUBGROUP}_multi_tissue_gene_expr_raw
 echo "-------------------------------------------------------------------------------"
 echo " "
done #for: SUBGROUP

head -10 $output_raw_dir/ever_tobacco_user_birthyear_groupingBIN4_multi_tissue_gene_expr_raw.cell_type_results.txt


