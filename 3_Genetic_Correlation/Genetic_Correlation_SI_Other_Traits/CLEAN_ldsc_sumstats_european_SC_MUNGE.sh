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
PHENOTYPE="SC"
GWAS="Liu_et_al_2019"

#Set directories
echo "Define directories."
root_dir="/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc"
workspace_dir="${root_dir}/Other_Sumstats/${PHENOTYPE}"
echo "Done."

#Settings (Cont')
echo "Settings."
REF_ALLELES="$root_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr/w_hm3.snplist"


echo ""
head -5 $workspace_dir/${GWAS}/SmokingCessation.txt
echo ""
echo "Count number of lines: "
wc -l $workspace_dir/${GWAS}/SmokingCessation.txt

munge_sumstats.py \
        --sumstats $workspace_dir/${GWAS}/SmokingCessation.txt \
        --chunksize 500000 \
        --snp RSID  \
        --a1 ALT \
        --a2 REF \
        --p PVALUE \
        --N-col N \
        --maf-min ${MAF_thr} \
        --signed-sumstats BETA,0 \
        --merge-alleles ${REF_ALLELES} \
        --out $workspace_dir/${GWAS}/munge_sumstats/${PHENOTYPE}_${GWAS}

echo " "
echo "Count number of remaining snps"
zcat $workspace_dir/${GWAS}/munge_sumstats/${PHENOTYPE}_${GWAS}.sumstats.gz | wc -l



