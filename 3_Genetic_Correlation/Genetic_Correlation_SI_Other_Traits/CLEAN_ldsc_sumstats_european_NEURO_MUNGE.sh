#Load the necessary modules (adjust based on your HPC system)
module load biocontainers/default
echo "Loading [biocontainers/default]: Done."
module load ldsc
echo "loading [ldsc]: done."

echo " "
echo "< Check whether [ldsc] is well loaded.>"
ldsc.py -h
echo " "
echo "< Check whether [munge_sumstats] is well loaded.>"
munge_sumstats.py -h
echo "Done."

echo " "
#Settings
echo "Settings."
ANCESTRY="european"
MAF_thr=0.01 # MAF threshold. 
PHENOTYPE="NEURO"
GWAS="Becker_et_al_2021"

#Set directories
echo "Define directories."
root_dir="/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc"
workspace_dir="${root_dir}/Other_Sumstats/${PHENOTYPE}"
echo "Done."

#Settings (Cont')
echo "Settings."
REF_ALLELES="$root_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr/w_hm3.snplist"


echo ""
head -5 $workspace_dir/${GWAS}/${PHENOTYPE}/${PHENOTYPE}1_excl_23andMe_single_gwide_sumstats.txt
echo ""
echo "Count number of lines: "
wc -l $workspace_dir/${GWAS}/${PHENOTYPE}/${PHENOTYPE}1_excl_23andMe_single_gwide_sumstats.txt

munge_sumstats.py \
        --sumstats $workspace_dir/${GWAS}/${PHENOTYPE}/${PHENOTYPE}1_excl_23andMe_single_gwide_sumstats.txt \
        --chunksize 500000 \
        --snp SNPID  \
        --a1 EFFECT_ALLELE \
        --a2 OTHER_ALLELE \
        --p PVALUE \
        --N-col N \
        --maf-min ${MAF_thr} \
        --signed-sumstats BETA,0 \
        --merge-alleles ${REF_ALLELES} \
        --out $workspace_dir/${GWAS}/munge_sumstats/${PHENOTYPE}_${GWAS}

echo " "
echo "Count number of remaining snps"
zcat $workspace_dir/${GWAS}/munge_sumstats/${PHENOTYPE}_${GWAS}.sumstats.gz | wc -l


echo " "
#Settings
echo "Settings."
ANCESTRY="european"
MAF_thr=0.01 # MAF threshold. 
PHENOTYPE="NEURO"
phenotype="Neuro"
GWAS="Gupta_et_al_2024"

#Set directories
echo "Define directories."
root_dir="/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc"
workspace_dir="${root_dir}/Other_Sumstats/${PHENOTYPE}"
echo "Done."

#Settings (Cont')
echo "Settings."
REF_ALLELES="$root_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr/w_hm3.snplist"


echo ""
head -5 $workspace_dir/${GWAS}/${phenotype}_MVP_UKB_sumstat_file
echo ""
echo "Count number of lines: "
wc -l $workspace_dir/${GWAS}/${phenotype}_MVP_UKB_sumstat_file

#NOTE: as the summary stat is VCF format we need to extract necessary columns.
awk 'BEGIN{OFS="\t"} 
     NR==1{ 
     print "ref","alt","rsid","N","beta","se","pval"; 
     next} 
     {
      print $5, $4, $1, $11, $7, $8, $9}' \
    $workspace_dir/${GWAS}/${phenotype}_MVP_UKB_sumstat_file \
    > $workspace_dir/${GWAS}/${PHENOTYPE}.txt

echo ""
echo "[${PHENOTYPE}]"
head -5 $workspace_dir/${GWAS}/${PHENOTYPE}.txt
echo ""
echo "Count number of lines: "
wc -l $workspace_dir/${GWAS}/${PHENOTYPE}.txt

munge_sumstats.py \
        --sumstats $workspace_dir/${GWAS}/${PHENOTYPE}.txt \
        --chunksize 500000 \
        --snp rsid  \
        --a1 alt \
        --a2 ref \
        --p pval \
        --N-col N \
        --maf-min ${MAF_thr} \
        --signed-sumstats beta,0 \
        --merge-alleles ${REF_ALLELES} \
        --out $workspace_dir/${GWAS}/munge_sumstats/${PHENOTYPE}_${GWAS}

echo " "
echo "Count number of remaining snps"
zcat $workspace_dir/${GWAS}/munge_sumstats/${PHENOTYPE}_${GWAS}.sumstats.gz | wc -l



