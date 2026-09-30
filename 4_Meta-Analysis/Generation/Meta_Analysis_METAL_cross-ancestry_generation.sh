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
type="generation"
pair="allancestry"
SUBGROUP_LIST1=($(for i in {4..7}; do echo "birthyear_bin$i"; done))
SUBGROUP_LIST2=($(for i in {8..11}; do echo "birthyear_bin$i"; done))
SUBGROUP_LIST3=($(for i in {12..15}; do echo "birthyear_bin$i"; done))
SUBGROUP_LIST4=($(for i in {16..19}; do echo "birthyear_bin$i"; done))

# Ensure the output folder exists
#mkdir -p "${working_dir}/Meta_Analysis/${SUBGROUP}"

# Loop over BIN indices
for bin in {1..4}; do

# Choose the right subgroup list dynamically
LISTNAME="SUBGROUP_LIST${bin}"

  # Make a nameref "CUR_LIST" that points to the array named in $LISTNAME
  declare -n CUR_LIST="${LISTNAME}"

SCRIPT=${working_dir}/${mod}/METAL_${SUBGROUP}_${type}${bin}_${pair}_script.txt
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
PVALUE P-value

# -------------------------------------------------------------------------------------------
EOF
# -------------------------------------------------------------------------------------------


##NOTE: If the input file formats are the same across input files, you do not need to specify the file formats again.
# -------------------------------------------------------------------------------------------
# Append PROCESS commands to the script using a loop
for LIST in "${CUR_LIST[@]}"; do

  data_dir=${working_dir}/${mod}/METAL_${SUBGROUP}_${LIST}_${pair}_MAresults_cleaned.txt

  echo "PROCESS ${data_dir}" >> ${SCRIPT}

done # LIST
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Append remaining commands to the script
cat <<EOF >> ${SCRIPT}

# ==========================================================================================

# Optional: To save the results to a file with a specific name
#OUTFILE METAL_${ANCESTRY}_${SUBGROUP}_MAresults .txt
OUTFILE ${working_dir}/${mod}/METAL_${SUBGROUP}_${type}${bin}_${pair}_MAresults .txt

# Optional: To restrict the output to only markers that have at least a specific number of individuals analyzed (or weight),
# use a command like the following:
# MINWEIGHT 10000 


# Perform the meta-analysis and calculate heterogeneity statistics
ANALYZE HETEROGENEITY


# ================================= [END: Write a script] ==================================
EOF
# ------------------------------------------------------------------------------------------

done #bin

# Change to output directory
#cd ${working_dir}"/Meta_Analysis/"${SUBGROUP}
# ------------------------------------- RUNNING SCRIPT -------------------------------------
# Loop over BIN indices
for bin in {1..4}; do
 echo ""
 echo "===================="
 echo "Generation ${bin} running..."
 echo "===================="
 echo ""
 ${software_meta_dir}/metal ${working_dir}/${mod}/METAL_${SUBGROUP}_${type}${bin}_${pair}_script.txt
 echo "===================="
done #bin
# ------------------------------------------------------------------------------------------

for i in {1..4}; do

echo ""
echo "--------------------------------------[${type}${i}]---------------------------------------"
echo ""

head -10 ${working_dir}/${mod}/METAL_${SUBGROUP}_${type}${i}_${pair}_MAresults1.txt
wc -l ${working_dir}/${mod}/METAL_${SUBGROUP}_${type}${i}_${pair}_MAresults1.txt

echo ""
echo "---------------------------------------------------------------------------------------"
echo ""
done #i


