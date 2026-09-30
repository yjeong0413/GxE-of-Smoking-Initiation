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
PHENOTYPE="MDD"
phenotype="MDD"
GWAS="PGC"

#Set directories
echo "Define directories."
root_dir="/YOUR PATH HERE/"
workspace_dir="${root_dir}/Other_Sumstats/${PHENOTYPE}"
echo "Done."

#Settings (Cont')
echo "Settings."
REF_ALLELES="$root_dir/${ANCESTRY}/Source_Files/eur_w_ld_chr/w_hm3.snplist"


echo "[${GWAS}]"
head -5 $workspace_dir/${GWAS}/PGC_${phenotype}.tsv
echo ""

echo "[PGC_${phenotype}.tsv]"
grep -v "^##" $workspace_dir/${GWAS}/PGC_${phenotype}.tsv | head -5

echo ""
echo "Count number of lines: "
wc -l $workspace_dir/${GWAS}/PGC_${phenotype}.tsv

#NOTE: as the summary stat is VCF format we need to extract necessary columns.
awk 'BEGIN{OFS="\t"} 
     /^##/ {next} 
     /^#CHROM/ { 
     print "ref","alt","rsid","N","beta","se","pval"; 
     next} 
     {N = $13 + $14;
      print $5, $4, $3, N, $6, $7, $8}' \
    $workspace_dir/${GWAS}/PGC_${phenotype}.tsv \
    > $workspace_dir/${GWAS}/${phenotype}.txt

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



