# Choose ancestry
#ANCESTRY="european"

# Set directories
echo "Set directories."
root_dir="/YOUR PATH HERE/"
working_dir=${root_dir}"/Cross-Ancestry_Meta-Analysis"

software_meta_dir="/YOUR PATH HERE/metal/generic-metal"
echo "Done."

# Change directory to [working_dir]
#cd ${working_dir}
pwd

SUBGROUP="cross-ancestry"
mod="birthyear"
pair="allancestry"
#SUBGROUP_LIST=("female" "male")

# loop through all combinations (i < j)
for bin_1 in {1..3}; do
    for bin_2 in $(seq $((bin_1+1)) 4); do
    
SCRIPT=${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}_script.txt

# ------------------------------------- WRITING SCRIPT -------------------------------------
cat <<EOF > ${SCRIPT}
# ================================ [START: Write a script] =================================

# Classical approach, uses effect size estimates and standard errors
SCHEME STDERR # Perform a meta-analysis weighted by the inverse of the squared standard erros (i.e., inverse variance weighting) It utilizes the effect size estimates (e.g., beta coefficients) and their corresponding standard errors to determine the weight each study contributes to the overall meta-analysis results.
#NOTE:
## WEIGHT (sample size) is optional if you’re using SCHEME STDERR, but can help with heterogeneity checks or SCHEME SAMPLESIZE.
### (When SCHEME STDERR and WEIGHT are both specified, METAl prioritizes the SCHEME STDERR command for calculating the primary meta-analysis weights.)
## FREQ is optional in METAL (used for QC/annotation but not for weighting unless SCHEME is frequency-based — which is rare).


# Turn on below if genomic control is needed.
# GENOMICCONTROL ON  # Uncomment and potentially specify a value if desired (e.g., GENOMICCONTROL 1.05)

# To help identify allele flips and check consistency, it can be useful to track allele frequencies in the meta-analysis. 
# Turn on below if needed.
AVERAGEFREQ ON
MINMAXFREQ ON
FREQLABEL Freq1  # Specify the header for your allele frequency column if it's not 'Freq1' or similar by default.

#summarize study specific results
#VERBOSE ON

# =============================== DESCRIBE THE INPUT FILE ===================================
## First input file format
MARKER MarkerName
ALLELE Allele1 Allele2 # effect allele (Allele1) and other_allele (Allele2)
FREQ Freq1  # If you have an allele frequency column (for the effect allele) in your input files, specify it here.
WEIGHT N_tot #Using N_tot as weight. If SCHEME STDERR is specified above, SCHEME STDERR will be prioritized. However, [N_tot] will still be used to calculate other statistics (e.g., Q and I^2).
EFFECT Effect
STDERR StdErr 
PVALUE P-Value

# -------------------------------------------------------------------------------------------
EOF
# -------------------------------------------------------------------------------------------


##NOTE: If the input file formats are the same across input files, you do not need to specify the file formats again.
# -------------------------------------------------------------------------------------------
# Append PROCESS commands to the script using a loop
  data_dir=${working_dir}/${mod}
  
  echo "PROCESS ${data_dir}/METAL_${SUBGROUP}_generation${bin_1}_${pair}_MAresults_cleaned_temp.txt" >> ${SCRIPT}
  echo "PROCESS ${data_dir}/METAL_${SUBGROUP}_generation${bin_2}_${pair}_MAresults_cleaned_temp.txt" >> ${SCRIPT}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Append remaining commands to the script
cat <<EOF >> ${SCRIPT}

# ==========================================================================================

# Optional: To save the results to a file with a specific name
#OUTFILE METAL_${ANCESTRY}_${SUBGROUP}_MAresults .txt
OUTFILE ${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}_MAresults .txt

# Optional: To restrict the output to only markers that have at least a specific number of individuals analyzed (or weight),
# use a command like the following:
# MINWEIGHT 10000 


# Perform the meta-analysis and calculate heterogeneity statistics
ANALYZE HETEROGENEITY


# ================================= [END: Write a script] ==================================
EOF
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

    echo "Created: ${SCRIPT}"
    echo ""
    
    done # bin_2
done # bin_1

# ------------------------------------- RUNNING SCRIPT -------------------------------------
for bin_1 in {1..3}; do
  for bin_2 in $(seq $((bin_1+1)) 4); do
    script_path=${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}_script.txt

    echo "=============================================="
    echo "Running METAL for BIN${bin_1} & BIN${bin_2}"
    echo "Script: $script_path"
    echo "----------------------------------------------"

    "${software_meta_dir}/metal" "$script_path"
    
    echo ""
    echo "Finished BIN${bin_1}${bin_2}"
    echo ""
  done # bin_2
done #bin_1
    

for bin_1 in {1..3}; do
  for bin_2 in $(seq $((bin_1+1)) 4); do
    output_path=${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}_MAresults1.txt

    echo "=============================================="
    echo "Checking METAL results for BIN${bin_1} & BIN${bin_2}"
    echo "Results: $output_path"
    echo "----------------------------------------------"

    head -10 ${output_path}    
    wc -l ${output_path}    
        
    echo ""
    echo "Finished BIN${bin_1}${bin_2}"
    echo ""
  done # bin_2
done #bin_1
    

bin_1="1"
bin_2="2"
bin_3="3"
bin_4="4"

SCRIPT=${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}${bin_3}${bin_4}_script.txt

# ------------------------------------- WRITING SCRIPT -------------------------------------
cat <<EOF > ${SCRIPT}
# ================================ [START: Write a script] =================================

# Classical approach, uses effect size estimates and standard errors
SCHEME STDERR # Perform a meta-analysis weighted by the inverse of the squared standard erros (i.e., inverse variance weighting) It utilizes the effect size estimates (e.g., beta coefficients) and their corresponding standard errors to determine the weight each study contributes to the overall meta-analysis results.
#NOTE:
## WEIGHT (sample size) is optional if you’re using SCHEME STDERR, but can help with heterogeneity checks or SCHEME SAMPLESIZE.
### (When SCHEME STDERR and WEIGHT are both specified, METAl prioritizes the SCHEME STDERR command for calculating the primary meta-analysis weights.)
## FREQ is optional in METAL (used for QC/annotation but not for weighting unless SCHEME is frequency-based — which is rare).


# Turn on below if genomic control is needed.
# GENOMICCONTROL ON  # Uncomment and potentially specify a value if desired (e.g., GENOMICCONTROL 1.05)

# To help identify allele flips and check consistency, it can be useful to track allele frequencies in the meta-analysis. 
# Turn on below if needed.
AVERAGEFREQ ON
MINMAXFREQ ON
FREQLABEL Freq1  # Specify the header for your allele frequency column if it's not 'Freq1' or similar by default.

#summarize study specific results
#VERBOSE ON

# =============================== DESCRIBE THE INPUT FILE ===================================
## First input file format
MARKER MarkerName
ALLELE Allele1 Allele2 # effect allele (Allele1) and other_allele (Allele2)
FREQ Freq1  # If you have an allele frequency column (for the effect allele) in your input files, specify it here.
WEIGHT N_tot #Using N_tot as weight. If SCHEME STDERR is specified above, SCHEME STDERR will be prioritized. However, [N_tot] will still be used to calculate other statistics (e.g., Q and I^2).
EFFECT Effect
STDERR StdErr 
PVALUE P-Value

# -------------------------------------------------------------------------------------------
EOF
# -------------------------------------------------------------------------------------------


##NOTE: If the input file formats are the same across input files, you do not need to specify the file formats again.
# -------------------------------------------------------------------------------------------
# Append PROCESS commands to the script using a loop
  data_dir=${working_dir}/${mod}
  
  echo "PROCESS ${data_dir}/METAL_${SUBGROUP}_generation${bin_1}_${pair}_MAresults_cleaned_temp.txt" >> ${SCRIPT}
  echo "PROCESS ${data_dir}/METAL_${SUBGROUP}_generation${bin_2}_${pair}_MAresults_cleaned_temp.txt" >> ${SCRIPT}
  echo "PROCESS ${data_dir}/METAL_${SUBGROUP}_generation${bin_3}_${pair}_MAresults_cleaned_temp.txt" >> ${SCRIPT}
  echo "PROCESS ${data_dir}/METAL_${SUBGROUP}_generation${bin_4}_${pair}_MAresults_cleaned_temp.txt" >> ${SCRIPT}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Append remaining commands to the script
cat <<EOF >> ${SCRIPT}

# ==========================================================================================

# Optional: To save the results to a file with a specific name
#OUTFILE METAL_${ANCESTRY}_${SUBGROUP}_MAresults .txt
OUTFILE ${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}${bin_3}${bin_4}_MAresults .txt

# Optional: To restrict the output to only markers that have at least a specific number of individuals analyzed (or weight),
# use a command like the following:
# MINWEIGHT 10000 


# Perform the meta-analysis and calculate heterogeneity statistics
ANALYZE HETEROGENEITY


# ================================= [END: Write a script] ==================================
EOF
# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------

    echo "Created: ${SCRIPT}"
    echo ""


# ------------------------------------- RUNNING SCRIPT -------------------------------------
    script_path=${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}${bin_3}${bin_4}_script.txt

    echo "=============================================="
    echo "Running METAL for BIN${bin_1} & BIN${bin_2} & BIN${bin_3} & BIN${bin_4}"
    echo "Script: $script_path"
    echo "----------------------------------------------"

    "${software_meta_dir}/metal" "$script_path"
    
    echo ""
    echo "Finished BIN${bin_1}${bin_2}${bin_3}${bin_4}"
    echo ""
    

head -10 ${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}${bin_3}${bin_4}_MAresults1.txt
echo
wc -l ${working_dir}/${mod}/METAL_${SUBGROUP}MA_generation${bin_1}${bin_2}${bin_3}${bin_4}_MAresults1.txt


