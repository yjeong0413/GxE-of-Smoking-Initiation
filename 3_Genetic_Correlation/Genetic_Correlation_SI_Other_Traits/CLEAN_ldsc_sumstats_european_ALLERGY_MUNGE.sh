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
PHENOTYPE="ALLERGY"
phenotype="ALLERGY"
GWAS="NEALE_v2"

#Set directories
echo "Define directories."
root_dir="/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc"
workspace_dir="${root_dir}/Other_Sumstats/${PHENOTYPE}"
echo "Done."

#Settings (Cont')
echo "Settings."
REF_ALLELES="$root_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr/w_hm3.snplist"


echo "[variants]"
head -5 ${root_dir}/Other_Sumstats//NEALE_v2_annotation/variants.tsv
echo ""
echo "[${GWAS}]"
head -5 $workspace_dir/${GWAS}/${phenotype}.gwas.imputed_v3.both_sexes.tsv
echo ""
echo "Count number of lines: "
wc -l $workspace_dir/${GWAS}/${phenotype}.gwas.imputed_v3.both_sexes.tsv

#NOTE: The order of variants in this file matches the order of variants in the results files described below. 
#      To join these annotations with a results file, either match on the "variant" field or simply paste the 
#      columns together (e.g. "paste variants.tsv K50.gwas.imputed_v3.both_sexes.tsv").


# combine two files
paste \
    <(cut -f4,5,6 ${root_dir}/Other_Sumstats//NEALE_v2_annotation/variants.tsv) \
    <(cut -f6,9,10,12 $workspace_dir/${GWAS}/${phenotype}.gwas.imputed_v3.both_sexes.tsv) \
    | sed '1s/n_complete_samples/N/' \
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



