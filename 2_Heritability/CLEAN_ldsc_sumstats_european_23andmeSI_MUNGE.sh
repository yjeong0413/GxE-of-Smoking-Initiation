#Load the necessary modules (adjust based on your HPC system)
module load biocontainers/default
echo "Loading [biocontainers/default]: Done.}"
module load ldsc
echo "loading [ldsc]: done."

echo " "
echo "< Check whether [ldsc] is well loaded.>"
ldsc.py -h
echo " "
echo "< Check whether [munge_sumstats] is well loaded.>"
munge_sumstats.py -h
echo "Done."

#Settings
echo "Settings."
ANCESTRY="european"
MAF_thr=0.01 # MAF threshold. 

echo " "
#Set directories
echo "Define directories."
workspace_dir="/YOUR PATH HERE/"
root_dir="${workspace_dir}/${ANCESTRY}/CLEAN"
dest_dir="${workspace_dir}/${ANCESTRY}/CLEAN/GWAS_QC"

#Settings
echo ""
echo "Settings."
REF_ALLELES="$workspace_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr/w_hm3.snplist"


echo " "
echo "Done."

# Check the data
echo "Fist 5 rows:"
head -5 $dest_dir/ever_tobacco_user/merged_CLEAN3.txt
echo " "
echo "Last 5 rows:"
tail -5 $dest_dir/ever_tobacco_user/merged_CLEAN3.txt
echo " "
echo "Total number of rows:"
wc -l $dest_dir/ever_tobacco_user/merged_CLEAN3.txt

echo " "
# FULL sample
subsample=("full")

for key in "${subsample[@]}"; do
echo "#----------------------------------------- [${key}] ------------------------------------------#"
echo " "
munge_sumstats.py \
        --sumstats $dest_dir/ever_tobacco_user/merged_CLEAN3.txt \
        --chunksize 500000 \
        --snp snpid \
        --a1 B \
        --a2 A \
        --p pval \
        --maf-min ${MAF_thr} \
        --signed-sumstats effect,0 \
        --N-cas-col N_cas \
        --N-con-col N_con \
        --merge-alleles ${REF_ALLELES} \
        --out $dest_dir/ever_tobacco_user/ever_tobacco_user

echo " "
echo "Count number of remaining snps"
zcat $dest_dir/ever_tobacco_user/ever_tobacco_user.sumstats.gz | wc -l

done #for: key

echo " "
#GENDER
subsample=("female" "male")

for key in "${subsample[@]}"; do
echo "#----------------------------------------- [${key}] ------------------------------------------#"
echo " "
echo " "
munge_sumstats.py \
        --sumstats $dest_dir/ever_tobacco_user_${key}/merged_CLEAN3.txt \
        --chunksize 500000 \
        --snp snpid \
        --a1 B \
        --a2 A \
        --p pval \
        --maf-min ${MAF_thr} \
        --signed-sumstats effect,0 \
        --N-cas-col N_cas \
        --N-con-col N_con \
        --merge-alleles ${REF_ALLELES} \
        --out $dest_dir/ever_tobacco_user_${key}/ever_tobacco_user_${key}


echo " "
echo "Count number of remaining snps"
zcat $dest_dir/ever_tobacco_user_${key}/ever_tobacco_user_${key}.sumstats.gz | wc -l

done #for: key

echo " "
#REGION
subsample=("MW" "NE" "SE" "SW" "W")

for key in "${subsample[@]}"; do
echo "#----------------------------------------- [${key}] ------------------------------------------#"
echo " "
echo " "
munge_sumstats.py \
        --sumstats $dest_dir/ever_tobacco_user_${key}/merged_CLEAN3.txt \
        --chunksize 500000 \
        --snp snpid \
        --a1 B \
        --a2 A \
        --p pval \
        --maf-min ${MAF_thr} \
        --signed-sumstats effect,0 \
        --N-cas-col N_cas \
        --N-con-col N_con \
        --merge-alleles ${REF_ALLELES} \
        --out $dest_dir/ever_tobacco_user_${key}/ever_tobacco_user_${key}

echo " "
echo "Count number of remaining snps"
zcat $dest_dir/ever_tobacco_user_${key}/ever_tobacco_user_${key}.sumstats.gz | wc -l

done #for: key

echo " "
#BIRTHYEAR
subsample=(birthyear_bin{4..19})

for key in "${subsample[@]}"; do
echo "#----------------------------------------- [${key}] ------------------------------------------#"
echo " "
echo " "
munge_sumstats.py \
        --sumstats $dest_dir/ever_tobacco_user_${key}/merged_CLEAN3.txt \
        --chunksize 500000 \
        --snp snpid \
        --a1 B \
        --a2 A \
        --p pval \
        --maf-min ${MAF_thr} \
        --signed-sumstats effect,0 \
        --N-cas-col N_cas \
        --N-con-col N_con \
        --merge-alleles ${REF_ALLELES} \
        --out $dest_dir/ever_tobacco_user_${key}/ever_tobacco_user_${key}


echo " "
echo "Count number of remaining snps"
zcat $dest_dir/ever_tobacco_user_${key}/ever_tobacco_user_${key}.sumstats.gz | wc -l

done #for: key

head -5 $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_groupingBIN1_cleaned_munge.txt

echo " "
#BIRTHYEAR
subsample=(groupingBIN{1..4})

for key in "${subsample[@]}"; do
echo "#----------------------------------------- [${key}] ------------------------------------------#"
echo " "
echo " "
munge_sumstats.py \
        --sumstats $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${key}_cleaned_munge.txt \
        --chunksize 500000 \
        --snp snpid \
        --a1 B \
        --a2 A \
        --p pval \
        --maf-min ${MAF_thr} \
        --signed-sumstats effect,0 \
        --N-cas-col N_cas \
        --N-con-col N_con \
        --merge-alleles ${REF_ALLELES} \
        --out $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${key}


echo " "
echo "Count number of remaining snps"
zcat $dest_dir/Meta_Analysis/birthyear/METAL_${ANCESTRY}_birthyear_MAresults_${key}.sumstats.gz | wc -l

done #for: key


