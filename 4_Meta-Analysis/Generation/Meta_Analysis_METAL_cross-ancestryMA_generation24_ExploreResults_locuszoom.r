R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#ANCESTRY="european"

cat("Set directories\n")
start_dir="/YOUR PATH HERE/"
#sumstat_dir=paste0(start_dir,"/",ANCESTRY)

root_dir=paste0(start_dir,"/WORKPLACE/ldsc/Cross-Ancestry_Meta-Analysis")
#dest_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC")

##MA_dir=paste0(dest_dir,"/Meta_Analysis")

cat("\n-------------------------------------------------------------------------------\n")

# check which packages are loaded
search()
#getAnywhere(select)

cat("\nLoad Libraries\n")
cat("\n-----------------------------------------------------------------------------------\n")
cat("*Load [data.table]\n")
library(data.table)
cat("\n")
cat("*Load [dplyr]\n")
library(dplyr)
cat("\n")
#cat("*Load [rlang]\n")
#library(rlang)
cat("\n")
cat("*Load [purrr]\n")
library(purrr)
cat("\n")
cat("*Load [tidyr]\n")
library(tidyr)
cat("\n")
cat("*Load [purr]\n")
library(purrr)
cat("\n")
cat("*Load [stringr]\n")
library(stringr)
cat("\n")
cat("*Load [ggplot2]\n")
library(ggplot2)
cat("\n")
cat("*Load [qqman]\n")
#install.packages("qqman")
library(qqman)
cat("\n")
cat("*Load [scales]\n")
library(scales)
cat("\n")
cat("*Load [patchwork]\n")
library(patchwork)
cat("\n")
cat("*Load [tidyverse]\n")
library(tidyverse)
#cat("\n")
#cat("\nLoad Libraries for Locus Zoom plot\n")
#cat("\n-----------------------------------------------------------------------------------\n")
#cat("\n")
#cat("*Load [locuszoomr]\n") #https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html
#library(locuszoomr)
#cat("\n")
#cat("*Load [EnsDb.Hsapiens.v86]\n")
#library(EnsDb.Hsapiens.v86) # For GRCh37, use [EnsDb.Hsapiens.v75]. # For GRCh38, use [].
cat("\n-----------------------------------------------------------------------------------\n")
search()

# Define subgroup
type <- "birthyear"
 ##type2: defined below
type3 <- "Generation"
category <- "generation"
subgroups <- paste0("generation",1:4)
subgroups

for (SUBGROUP in subgroups) {
    
    cat("\n")
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))

    file_path <- paste0(root_dir,"/",
                        type,"/METAL_cross-ancestry_",SUBGROUP,"_allancestry_MAresults_cleaned.txt")
                               
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("1) Read the data.\n")
    
    cleaned_data <- fread(file_path, sep = "\t", header = T)
    
    #Check [chr] and convert X,Y,M to 23, 24, 25 if exists
    cat("\n [cleaned_data$chr]: Before change \n")
    print(table(cleaned_data$chr))
    cleaned_data$chr <- ifelse(cleaned_data$chr=="X", 23, 
                        ifelse(cleaned_data$chr=="Y", 24,
                        ifelse(cleaned_data$chr %in% c("M"," MT"), 25, cleaned_data$chr)))
    
    cleaned_data$chr <- as.numeric(cleaned_data$chr) # ensure it's numeric
    cat("\n [cleaned_data$chr]: After change \n")
    print(table(cleaned_data$chr))
   
    # Save results as dynamic variable
    assign(paste0("cleaned_", SUBGROUP), cleaned_data, envir = .GlobalEnv)
    
    data_name <- paste0("cleaned_", SUBGROUP)
    
    cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

   cat("\n------------------------------------------------------------------------------------------\n")
   cat("\n")

} #for: SUBGROUP

type2="generation24"

# Read MA results
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
   data <- fread(file_path, sep = "\t", header = T)

#Check [chr] and convert X,Y,M to 23, 24, 25 if exists
cat("\n [data$chr]: Before change \n")
print(table(data$chr))
data$chr <- ifelse(data$chr=="X", 23, 
            ifelse(data$chr=="Y", 24,
            ifelse(data$chr %in% c("M"," MT"), 25, data$chr)))
    
data$chr <- as.numeric(data$chr) # ensure it's numeric
cat("\n [data$chr]: After change \n")
print(table(data$chr))

# assign data to [MA_results
data_name <- paste0("MA_results")
assign(data_name,data, envir = .GlobalEnv)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))


MA_results_sigHetPVal <- MA_results %>%
                          dplyr::filter(-log10(HetPVal) > -log10(5*10^{-8})) 
    cat(paste0("\nReturn first 5 rows of [MA_results_sigHetPVal]\n"))
    print(head(MA_results_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [MA_results_sigHetPVal]\n"))
    print(tail(MA_results_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(MA_results_sigHetPVal))
    cat("\nTabulate [chr]\n")
    print(table(MA_results_sigHetPVal$chr))
    cat("\n")

# read hg19 position 
file_path <- paste0(root_dir, "/",type,"/FUMA/MA_results_",type2,"_hg19_pos.txt")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
data_info_tempp2 <- fread(file_path, sep = "\t", header = F)

cat("- Assign column names.\n")
colnames(data_info_tempp2) <- c("chr", "start", "end", "snpid")

data_info_temp2 <- data_info_tempp2 %>% select(snpid, end)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(data_info_temp2, n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(data_info_temp2, n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(data_info_temp2))


# join [data_info_temp1] and [data_info_temp2] to generate [data_info]: hg19 position
data_info <- left_join(MA_results %>% select(chr, MarkerName), data_info_temp2, 
                       by=c("MarkerName"="snpid")) %>%
             select(chr, MarkerName, end) %>% rename(snpid=MarkerName , position=end)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(data_info, n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(data_info, n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(data_info))


# Read high-PIP variants
file_path <- paste0(root_dir, "/",type,
            "/sex_generation_moderated_loci_with_highPIP_variant_ids_rsid_filled_summary_with_locus_coords.tsv")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
   data_temp <- fread(file_path, sep = "\t", header = T)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(data_temp, n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(data_temp, n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(data_temp))

   data_fine_all <- data_temp %>%
                filter(moderated_domain == paste0("generation_moderated")) %>%
                select(moderated_domain, chr_hg19, highPIP_rsids)

    # assign data to [data_meta_loci]
    data_name <- paste0("data_high_pip_all")
    assign(data_name,data_fine_all, envir = .GlobalEnv)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

    # generate only one column containing the unique rsIDs
    data_fine_rsid <- data_fine_all %>%
                    separate_rows(highPIP_rsids, sep = ";") %>%
                    transmute(MarkerName = trimws(highPIP_rsids)) %>%
                    filter(!is.na(MarkerName), MarkerName != "") %>%
                    distinct()

    # assign data to [data_meta_loci]
    data_name <- paste0("data_high_pip_rsid")
    assign(data_name,data_fine_rsid, envir = .GlobalEnv)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))



# Read MA results
file_path <- paste0(root_dir, "/",type,"/meta_locus_table_for_locuszoom_with_original_hg19.tsv")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
   data_temp <- fread(file_path, sep = "\t", header = T)

   data <- data_temp %>%
                        filter(locus_class == paste0(category,"_moderated"),
                              .data[[paste0(category,"_contrast")]] == "G2_vs_G4")

    # assign data to [data_meta_loci]
    data_name <- paste0("data_meta_loci")
    assign(data_name,data, envir = .GlobalEnv)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

FUMA_job="FUMA_job737777_bin24_05-15-2026"

# Read MA results
file_path <- paste0(root_dir, "/",type,"/FUMA/Results/",FUMA_job,"/GenomicRiskLoci.txt")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
   data_temp <- fread(file_path, sep = "\t", header = T)

    data <- data_temp %>%
                select(chr, rsID, pos, start, end, GenomicLocus)

    cat("- Assign column names.\n")
    colnames(data) <- c("chr", "leadSNPs", "pos_hg19", "start", "end", "GenomicRiskLoci")

    # merge with [data_meta_loci]
    data_merged <- left_join(data, 
                             data_meta_loci %>% select(meta_start_hg19, meta_end_hg19,
                                                       original_chr_hg19,
                                                       original_start_hg19, original_end_hg19),
                             by = c("chr"="original_chr_hg19",
                                    "start"="original_start_hg19",
                                    "end"="original_end_hg19"))

    # assign data to [FUMA_results_GenomicRiskLoci]
    data_name <- paste0("FUMA_results_GenomicRiskLoci")
    assign(data_name,data_merged, envir = .GlobalEnv)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

ls()

# Read MA results
file_path <- paste0(root_dir, "/",type,"/FUMA/Results/",FUMA_job,"/genes.txt")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
   data_temp <- fread(file_path, sep = "\t", header = T)

    data <- data_temp %>%
                select(chr, symbol, type, start, end)

    cat("- Assign column names.\n")
    colnames(data) <- c("chr", "symbol", "type", "start_gene", "end_gene")

    # assign data to [MA_results
    data_name <- paste0("FUMA_results")
    assign(data_name,data, envir = .GlobalEnv)

    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

ls()

############################################################
### LocusZoom-style plots by FUMA Genomic Risk Locus
### Female vs Male example
############################################################
############################################################
### Settings
############################################################

# Group
group_1 <- "generation2"
group_2 <- "generation4"

# Legend Names
Group_1 <- "Baby Boomers"
Group_2 <- "Millennials"

# Genome-wide significance threshold
genome_sig <- 5e-8

# Extra base-pair window around each genomic risk locus
bp_extra <- 0 #10000

# Container for snp with lowest HetPVal in each locus
lowest_het_snp_list <- list()
# Container for snp with moderated loci
sig_snps_locus_list <- list()

# Chromosome
for (CHR in c(6:7,10:11,14,16)) {

cat("\n*--------------------------[Chromosome ",CHR,"]---------------------------*\n")
    

############################################################
### Assign MA data to consider
############################################################

data_temp1 <- get(paste0("cleaned_", group_1))
print(dim(data_temp1))

data_temp2 <- get(paste0("cleaned_", group_2))
print(dim(data_temp2))

data_MA <- get(paste0("MA_results"))
print(dim(data_MA))


############################################################
### Filter two MA results to have only common SNPs with MA results
############################################################

data_1 <- data_temp1[data_temp1$MarkerName %in% data_MA$MarkerName, ]
print(dim(data_1))

data_2 <- data_temp2[data_temp2$MarkerName %in% data_MA$MarkerName, ]
print(dim(data_2))

cat("\n[data_1] \n")
print(head(data_1, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(data_1))

cat("\n[data_2] \n")
print(head(data_2, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(data_2))


############################################################
### Extract effect and p-value from the two subgroup results
############################################################

new_name_1  <- paste0("Effect_", group_1)
new_name_1p <- paste0("Pvalue_", group_1)

new_name_2  <- paste0("Effect_", group_2)
new_name_2p <- paste0("Pvalue_", group_2)

merged_MA_temp <- data_MA %>%
  dplyr::select(MarkerName, HetPVal) %>%
  left_join(
    dplyr::select(data_1, MarkerName, Effect_new, "P-value"), by = c("MarkerName")
  ) %>%
  dplyr::rename(
    !!new_name_1 := Effect_new,
    !!new_name_1p := "P-value"
  ) %>%
  left_join(
    dplyr::select(data_2, MarkerName, Effect_new, "P-value"), by = c("MarkerName")
  ) %>%
  dplyr::rename(
    !!new_name_2 := Effect_new,
    !!new_name_2p := "P-value"
  )

cat("\n[merged_MA_temp]\n")
print(head(merged_MA_temp, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA_temp))


############################################################
### Merge chromosome and position (hg19) information, then filter by CHR
############################################################

merged_MA <- merged_MA_temp %>%
  left_join(data_info, by = c("MarkerName" = "snpid")) %>%
  dplyr::filter(chr == CHR)

cat("\n[merged_MA]\n")
print(head(merged_MA, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA))


############################################################
### Identify SNPs with at least one subgroup p-value < 5e-8
############################################################

merged_MA_sigEffect <- merged_MA %>%
  dplyr::filter(
    (!!sym(new_name_1p)) < genome_sig |
    (!!sym(new_name_2p)) < genome_sig
  )

cat("\nReturn first 5 rows of [merged_MA_sigEffect]\n")
print(head(merged_MA_sigEffect, n = 5))

cat("\nReturn last 5 rows of [merged_MA_sigEffect]\n")
print(tail(merged_MA_sigEffect, n = 5))

cat("\nSize of [merged_MA_sigEffect]: (# of rows, # of columns)\n")
print(dim(merged_MA_sigEffect))
cat("\n")


############################################################
### Identify SNPs with HetPVal < 5e-8
############################################################

merged_MA_sigHetPVal <- merged_MA %>%
  dplyr::filter(HetPVal < genome_sig)

cat("\nReturn first 5 rows of [merged_MA_sigHetPVal]\n")
print(head(merged_MA_sigHetPVal, n = 5))

cat("\nReturn last 5 rows of [merged_MA_sigHetPVal]\n")
print(tail(merged_MA_sigHetPVal, n = 5))

cat("\nSize of [merged_MA_sigHetPVal]: (# of rows, # of columns)\n")
print(dim(merged_MA_sigHetPVal))
cat("\n")


############################################################
### Get intersection of [significant effect] and [significant HetPVal SNPs]
############################################################

merged_MA_sigEffectHetPVal <- merged_MA_sigEffect %>%
  dplyr::filter(MarkerName %in% merged_MA_sigHetPVal$MarkerName)

cat("\nReturn first 5 rows of [merged_MA_sigEffectHetPVal]\n")
print(head(merged_MA_sigEffectHetPVal, n = 5))

cat("\nReturn last 5 rows of [merged_MA_sigEffectHetPVal]\n")
print(tail(merged_MA_sigEffectHetPVal, n = 5))

cat("\nSize of [merged_MA_sigEffectHetPVal]: (# of rows, # of columns)\n")
print(dim(merged_MA_sigEffectHetPVal))

cat("\n[chr]\n")
print(table(merged_MA_sigEffectHetPVal$chr))


############################################################
### Check how many significant SNPs fall into each FUMA Genomic Risk Locus
############################################################

genomicriskloci <- FUMA_results_GenomicRiskLoci %>% dplyr::filter(chr == CHR)

cat("\n[Identified Genomic Risk Loci]\n")
print(head(genomicriskloci, n = 5))

snps_to_check <- merged_MA_sigEffectHetPVal$position

genomicriskloci$SNP_count <- sapply(seq_len(nrow(genomicriskloci)), function(i) {
  sum(
    snps_to_check >= genomicriskloci$meta_start_hg19[i] &
      snps_to_check <= genomicriskloci$meta_end_hg19[i]
  )
})

cat("\n[Number of SNPs in merged_MA_sigEffectHetPVal in the identified Genomic Risk Loci]\n")
print(genomicriskloci)


############################################################
### Prepare data for effect-size flipping
############################################################

# Dynamic column names
col1 <- paste0("Effect_", group_1)
col2 <- paste0("Effect_", group_2)

# Use group 1 effect as reference direction
merged_MA$compare_effect <- merged_MA[[col1]]

# Flip signs so group 1 effect is always positive
merged_MA_flip <- merged_MA %>%
  mutate(
    !!sym(col1) := if_else(compare_effect < 0, -!!sym(col1), !!sym(col1)),
    !!sym(col2) := if_else(compare_effect < 0, -!!sym(col2), !!sym(col2))
  )

cat("\nReturn first 5 rows of [merged_MA]\n")
print(head(merged_MA, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA))

cat("\nReturn first 5 rows of [merged_MA_flip]\n")
print(head(merged_MA_flip, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA_flip))


############################################################
### Keep only genomic risk loci with at least one SNP to highlight
############################################################

genomicriskloci_to_plot <- genomicriskloci %>% dplyr::filter(SNP_count > 0)

cat("\n[Genomic Risk Loci to Plot]\n")
print(genomicriskloci_to_plot)

if (nrow(genomicriskloci_to_plot) == 0) {
  cat("No genomic risk loci contain SNPs from merged_MA_sigEffectHetPVal on CHR", CHR, 
      ". Skipping this chromosome.\n")
  next
}


############################################################
### Loop over each Genomic Risk Locus and create separate plots
############################################################

for (i in seq_len(nrow(genomicriskloci_to_plot))) {

  locus_i <- genomicriskloci_to_plot[i, ]

  ############################################################
  ### Define locus ID
  ############################################################

    locus_id <- i
    locus_start_pos <- locus_i$meta_start_hg19
    locus_end_pos <- locus_i$meta_end_hg19

  cat("\n====================================================\n")
  cat("Plotting CHR", CHR, "Genomic Risk Locus", locus_id, "\n")
  cat("Locus start:", locus_start_pos, "\n")
  cat("Locus end:", locus_end_pos, "\n")
  cat("====================================================\n")

  ############################################################
  ### Get hg38 and hg 19 position of FUMA lead SNP for plot title
  ############################################################
  lead_snp <- as.character(locus_i$leadSNPs)
  lead_pos_hg19 <- locus_i$pos_hg19
  lead_info <- data_temp1 %>%
                    dplyr::filter(MarkerName == lead_snp) %>%
                    dplyr::slice(1)
  lead_chr <- lead_info$chr
  lead_pos_hg38 <- lead_info$position
    
  cat("\n====================================================\n")
  cat("Lead SNP:", lead_snp, "\n")
  cat("Lead SNP chromosome:", lead_chr, "\n")
  cat("Lead SNP position (hg38):", lead_pos_hg38, "\n")
  cat("Lead SNP position (hg19):", lead_pos_hg19, "\n")
  cat("====================================================\n")


  ############################################################
  ### Define locus-specific plotting window
  ############################################################

  bp_min <- locus_i$meta_start_hg19 - bp_extra
  bp_max <- locus_i$meta_end_hg19 + bp_extra

  cat("\nPlotting window:\n")
  cat("bp_min =", bp_min, "\n")
  cat("bp_max =", bp_max, "\n")

  ############################################################
  ### Significant SNPs inside this specific locus
  ############################################################

  sig_snps_locus <- merged_MA_sigEffectHetPVal %>%
    dplyr::filter(
      position >= locus_i$meta_start_hg19,
      position <= locus_i$meta_end_hg19
    )
    
  cat("\n[sig_snps_locus]:\n")
  print(head(sig_snps_locus))
  
  # to save
  sig_snps_locus_save <- sig_snps_locus %>%
                            dplyr::select(MarkerName, HetPVal, chr, position) %>%
                            dplyr::mutate(moderation = type2,
                                          pos_hg19 = position)
  cat("\n[sig_snps_locus_save]:\n")
  print(head(sig_snps_locus_save))
    
  sig_snps_locus_list[[paste0("chr", CHR, "_locus", locus_id)]] <- sig_snps_locus_save
    
  # keep [MarkerName] only
  snp_mark_locus <- sig_snps_locus$MarkerName
    
  cat("\nNumber of highlighted SNPs in this locus:", length(snp_mark_locus), "\n")

  if (length(snp_mark_locus) == 0) {
    cat("Skipping this locus because no SNPs are available to highlight.\n")
    next
  }
    
############################################################
### High-PIP SNPs inside this specific locus
############################################################

high_pip_snps_locus <- merged_MA_flip %>%
  dplyr::filter(
    MarkerName %in% data_high_pip_rsid$MarkerName,
    position >= bp_min,
    position <= bp_max
  )

high_pip_mark_locus <- high_pip_snps_locus$MarkerName

cat("\nNumber of high-PIP SNPs in this locus:",
    length(high_pip_mark_locus), "\n")
    

  ############################################################
  ### SNP with lowest HetPVal inside this specific locus
  ############################################################

  # SNP with the lowest HetPVal in this locus
  lowest_het_snp_locus <- sig_snps_locus %>%
    dplyr::slice_min(
      order_by = HetPVal,
      n = 1,
      with_ties = FALSE
    )%>%
    dplyr::mutate(
      GenomicRiskLocus = locus_id
    )
    
  cat("\nSNP with lowest HetPVal in this locus:\n")
  print(lowest_het_snp_locus)

    
  # add start and end position to [lowest_het_snp_locus]
  lowest_het_snp_locus_new <- lowest_het_snp_locus %>%
                                    mutate(meta_start_hg19 = locus_start_pos,
                                          meta_end_hg19 = locus_end_pos)

  # Save to full list
  lowest_het_snp_list[[paste0("chr", CHR, "_locus", locus_id)]] <- lowest_het_snp_locus_new

  ############################################################
  ### Filter SNPs to locus window
  ############################################################

  merged_MA_small <- merged_MA_flip %>%
    dplyr::filter(
      position >= bp_min,
      position <= bp_max
    )

  cat("\nReturn first 5 rows of [merged_MA_small]\n")
  print(head(merged_MA_small, n = 5))

  cat("\nReturn last 5 rows of [merged_MA_small]\n")
  print(tail(merged_MA_small, n = 5))

  cat("\nSize of [merged_MA_small]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small))
  cat("\n")


  ############################################################
  ### Filter SNPs to highlight
  ############################################################

  merged_MA_small_mark <- merged_MA_flip %>%
    dplyr::filter(MarkerName %in% snp_mark_locus)

  cat("\nReturn first 5 rows of [merged_MA_small_mark]\n")
  print(head(merged_MA_small_mark, n = 5))

  cat("\nReturn last 5 rows of [merged_MA_small_mark]\n")
  print(tail(merged_MA_small_mark, n = 5))

  cat("\nSize of [merged_MA_small_mark]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small_mark))
  cat("\n")


  ############################################################
  ### Create long-form data frame for Panel A
  ############################################################

  merged_MA_small_long <- merged_MA_small %>%
    dplyr::select(MarkerName, chr, position, !!sym(col1), !!sym(col2)) %>%
    pivot_longer(
      cols = c(!!sym(col1), !!sym(col2)),
      names_to = "group",
      values_to = "Effect"
    )

  cat("\nReturn first 5 rows of [merged_MA_small_long]\n")
  print(head(merged_MA_small_long, n = 5))

  cat("\nSize of [merged_MA_small_long]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small_long))


  ############################################################
  ### Mark SNPs to highlight
  ############################################################

  merged_MA_small_mark_long <- merged_MA_small_long %>%
    dplyr::filter(MarkerName %in% snp_mark_locus)

  cat("\nReturn first 5 rows of [merged_MA_small_mark_long]\n")
  print(head(merged_MA_small_mark_long, n = 5))

  cat("\nSize of [merged_MA_small_mark_long]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small_mark_long))

  ############################################################
  ### Create high-PIP long data: Mark high-PIP SNPs
  ############################################################
  merged_MA_highPIP <- merged_MA_small %>%
    dplyr::filter(MarkerName %in% high_pip_mark_locus)

  merged_MA_highPIP_long <- merged_MA_small_long %>%
    dplyr::filter(MarkerName %in% high_pip_mark_locus)
    

  ############################################################
  ### Panel A: Marginal Genetic Effects Plot
  ############################################################

  effect_plot <- ggplot() +
    geom_segment(data = merged_MA_small,
                 aes(x = position,
                     xend = position,
                     y = !!sym(col1),
                     yend = !!sym(col2)),
                 linewidth = 0.25, alpha = 0.5, color = "gray75") +
    geom_point(data = dplyr::filter(merged_MA_small_long, group == col1),
               aes(position, Effect, color = paste0(Group_1, "_gray")),
               shape = 15, size = 2.3) +
    geom_point(data = dplyr::filter(merged_MA_small_long, group == col2),
               aes(position, Effect, color = paste0(Group_2, "_gray")),
               shape = 15, size = 2.3) +
    geom_segment(data = merged_MA_small_mark,
                 aes(x = position,
                     xend = position,
                     y = !!sym(col1),
                     yend = !!sym(col2)),
                 linewidth = 0.8,alpha = 0.5,color = "purple",show.legend = FALSE) +
    geom_hline(yintercept = 0,color = "gray",linewidth = 0.4,linetype = "dashed") +
    geom_point(data = dplyr::filter(merged_MA_small_mark_long, group %in% c(col1, col2)),
               aes(position, Effect, color = group),
               shape = 15,size = 2.2) +
    # High-PIP effect points
    geom_point(data = dplyr::filter(merged_MA_highPIP_long, group == col1),
               aes(position, Effect, color = "High-PIP: Group 1"),
               shape = 15, size = 2.8) +
    geom_point(data = dplyr::filter(merged_MA_highPIP_long, group == col2),
               aes(position, Effect, color = "High-PIP: Group 2"),
               shape = 15, size = 2.8) +    
    # Connect high-PIP SNP effects
    geom_segment(data = merged_MA_highPIP,
                 aes(x = position,
                     xend = position,
                     y = !!sym(col1),
                     yend = !!sym(col2)),
                 linewidth = 0.8, color = "#FFD600", show.legend = FALSE) +
    scale_color_manual(name = "",
                       values = c(
                       setNames("gray30", paste0(Group_1, "_gray")),
                       setNames("red", col1),
                       setNames("gray60", paste0(Group_2, "_gray")),
                       setNames("blue", col2),
                        "High-PIP: Group 1" = "#FFD600",
                        "High-PIP: Group 2" = "#FFD600"),
                      breaks = c(paste0(Group_1, "_gray"), col1,
                                 paste0(Group_2, "_gray"), col2,
                                 "High-PIP: Group 1"),
                      labels = c("",paste0(Group_1,"     "),
                                 "",paste0(Group_2,"     "),"High-PIP")) +
    guides(color = guide_legend(byrow = TRUE,nrow = 1,
                                keywidth = unit(0, "lines"))) +
    scale_x_continuous(labels = label_comma()) +
    coord_cartesian(xlim = c(bp_min, bp_max)) +
    labs(x = NULL, y = "Marginal Effects",
      title = paste0(Group_1, " vs ", Group_2, 
                     ", Chr", lead_chr, ":", locus_start_pos, "-", locus_end_pos)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 16),
      legend.spacing.x = unit(0, "cm"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.ticks.x = element_blank()
    )

  print(effect_plot)


  ############################################################
  ### Panel B: HetPVal Plot
  ############################################################

  df <- merged_MA_small %>%
    dplyr::mutate(
      HetPval = as.numeric(HetPVal),
      logHetP = -log10(HetPval)
    ) %>%
    dplyr::filter(!is.na(logHetP), is.finite(logHetP))

  hetpval_plot <- ggplot(df, aes(x = position, y = logHetP)) +
    geom_point(shape = 21,fill = "gray80",color = "black",stroke = 0.5,size = 2) +
    geom_hline(yintercept = -log10(genome_sig),
               linetype = "dashed",color = "gray60",linewidth = 0.4) +
    scale_x_continuous(labels = label_comma()) +
    coord_cartesian(xlim = c(bp_min, bp_max)) +
    labs(x = NULL,
         y = expression(-log[10]("HetP-value")),
         title = "") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.ticks.x = element_blank()
    ) +
    geom_point(data = df %>% dplyr::filter(MarkerName %in% snp_mark_locus),
      aes(position, logHetP),
      shape = 21,fill = "deepskyblue",color = "black",stroke = 0.5,size = 2
    )

  print(hetpval_plot)


  ############################################################
  ### Panel C: Mapped Genes
  ############################################################

  data_locus <- FUMA_results %>%
    dplyr::filter(chr == CHR) %>%
    dplyr::filter(type == "protein_coding") %>%
    dplyr::filter(
      end_gene >= bp_min,
      start_gene <= bp_max
    ) %>%
    mutate(
      start_mb = start_gene,
      end_mb = end_gene,
      plot_start = pmax(start_mb, bp_min),
      plot_end   = pmin(end_mb, bp_max),
      label_x = (plot_start + plot_end) / 2
    ) %>%
    arrange(start_mb) %>%
    mutate(y = row_number())

  data_name <- paste0("genes_chr", CHR, "_locus", locus_id)
  assign(data_name, data_locus)

  cat("\nReturn first 5 rows of mapped genes for this locus\n")
  print(head(get(data_name), n = 5))

  cat("\nReturn last 5 rows of mapped genes for this locus\n")
  print(tail(get(data_name), n = 5))

  cat("\nSize of mapped genes table: (# of rows, # of columns)\n")
  print(dim(get(data_name)))


  if (nrow(data_locus) > 0) {

    gene_track_chr <- ggplot(data_locus) +
      geom_segment(aes(x = start_mb,xend = end_mb,y = y,yend = y),
                   linewidth = 1.2,color = "black") +
      geom_text(aes(x = label_x,y = y,label = symbol),
                    size = 4,vjust = -0.6) +
      scale_y_continuous(breaks = NULL) +
      scale_x_continuous(labels = label_comma()) +
      coord_cartesian(
        xlim = c(bp_min, bp_max),
        ylim = c(min(data_locus$y) - 0.5, max(data_locus$y) + 0.5)) +
      labs(x = paste0("Base Pair Position on Chromosome ", CHR),
           y = "Mapped Genes") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 8)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)
      )

  } else {


    gene_track_chr <- ggplot() +
      coord_cartesian(xlim = c(bp_min, bp_max)) +
      scale_x_continuous(labels = label_comma()) +
      labs(x = paste0("Base Pair Position on Chromosome ", CHR),
           y = "Mapped Genes") +
      annotate("text",
                x = mean(c(bp_min, bp_max)),
                y = 1,
                label = "No protein-coding genes in this window",size = 5) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 8)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  print(gene_track_chr)


  ############################################################
  ### Combine Panels A, B, and C
  ############################################################

  combined_plot <- effect_plot / hetpval_plot / gene_track_chr +
    plot_layout(heights = c(1, 1, 1))

  print(combined_plot)


  ############################################################
  ### Save combined plot
  ############################################################

  output_file <- paste0(root_dir, "/", type,"/Output/Figures/locuszoom_chr", CHR,
                        "_locus", locus_id,"_", type2,"_sigHetPVal_FUMA.png")

  ggsave(output_file,combined_plot,width = 10,height = 9,dpi = 500)

  cat("\nSaved plot:\n")
  cat(output_file, "\n")
}
    
cat("\n*-----------------------------------------------------------------------------------*\n")

} # for: CHR


############################################################
### Combine lowest-HetPVal SNPs across all loci
############################################################

lowest_het_snp_all <- dplyr::bind_rows(
  lowest_het_snp_list,
  .id = "locus_list"
)

lowest_het_snp_all <- lowest_het_snp_all %>%
 dplyr::mutate(type=type2) %>% 
 dplyr::arrange(chr, position)

cat("\n[Lowest HetPVal SNP from each Genomic Risk Locus]\n")
print(lowest_het_snp_all)


# Save as txt
file_path = paste0(root_dir, "/", type,"/Output/", type3, "_Lowest_HetPval")
write.table(lowest_het_snp_all,
      file = paste0(file_path, "/LowestHetPValSNP_",type2,".txt"), 
      sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE
     )


############################################################
### Combine moderated variants across all loci
############################################################
sig_snps_locus_all <- dplyr::bind_rows(
  sig_snps_locus_list,
  .id = "locus_list"
)
cat("\n[Moderated SNP from each Meta Locus]: sig_snp_locus_all\n")
print(sig_snps_locus_all)

# Save as txt
# Define output directory
file_path <- paste0(root_dir, "/", type, "/Output/", type3, "_Moderated_Variants")

# Create directory if it does not exist
if (!dir.exists(file_path)) {dir.create(file_path, recursive = TRUE)}

# Save as txt
write.table(sig_snps_locus_all,
  file = paste0(file_path, "/Moderated_Variants_", type2, ".txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

############################################################
### LocusZoom-style plots by FUMA Genomic Risk Locus
### Female vs Male example
############################################################
############################################################
### Settings
############################################################

# Group
group_1 <- "generation2"
group_2 <- "generation4"

# Legend Names
Group_1 <- "Baby Boomers"
Group_2 <- "Millennials"

# Genome-wide significance threshold
genome_sig <- 5e-8

# Extra base-pair window around each genomic risk locus
bp_extra <- 50000

# Container for snp with lowest HetPVal in each locus
lowest_het_snp_list <- list()

# Chromosome
for (CHR in c(6:7,10:11,14,16)) {

cat("\n*--------------------------[Chromosome ",CHR,"]---------------------------*\n")
    

############################################################
### Assign MA data to consider
############################################################

data_temp1 <- get(paste0("cleaned_", group_1))
print(dim(data_temp1))

data_temp2 <- get(paste0("cleaned_", group_2))
print(dim(data_temp2))

data_MA <- get(paste0("MA_results"))
print(dim(data_MA))


############################################################
### Filter two MA results to have only common SNPs with MA results
############################################################

data_1 <- data_temp1[data_temp1$MarkerName %in% data_MA$MarkerName, ]
print(dim(data_1))

data_2 <- data_temp2[data_temp2$MarkerName %in% data_MA$MarkerName, ]
print(dim(data_2))

cat("\n[data_1] \n")
print(head(data_1, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(data_1))

cat("\n[data_2] \n")
print(head(data_2, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(data_2))


############################################################
### Extract effect and p-value from the two subgroup results
############################################################

new_name_1  <- paste0("Effect_", group_1)
new_name_1p <- paste0("Pvalue_", group_1)

new_name_2  <- paste0("Effect_", group_2)
new_name_2p <- paste0("Pvalue_", group_2)

merged_MA_temp <- data_MA %>%
  dplyr::select(MarkerName, HetPVal) %>%
  left_join(
    dplyr::select(data_1, MarkerName, Effect_new, "P-value"), by = c("MarkerName")
  ) %>%
  dplyr::rename(
    !!new_name_1 := Effect_new,
    !!new_name_1p := "P-value"
  ) %>%
  left_join(
    dplyr::select(data_2, MarkerName, Effect_new, "P-value"), by = c("MarkerName")
  ) %>%
  dplyr::rename(
    !!new_name_2 := Effect_new,
    !!new_name_2p := "P-value"
  )

cat("\n[merged_MA_temp]\n")
print(head(merged_MA_temp, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA_temp))


############################################################
### Merge chromosome and position (hg19) information, then filter by CHR
############################################################

merged_MA <- merged_MA_temp %>%
  left_join(data_info, by = c("MarkerName" = "snpid")) %>%
  dplyr::filter(chr == CHR)

cat("\n[merged_MA]\n")
print(head(merged_MA, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA))


############################################################
### Identify SNPs with at least one subgroup p-value < 5e-8
############################################################

merged_MA_sigEffect <- merged_MA %>%
  dplyr::filter(
    (!!sym(new_name_1p)) < genome_sig |
    (!!sym(new_name_2p)) < genome_sig
  )

cat("\nReturn first 5 rows of [merged_MA_sigEffect]\n")
print(head(merged_MA_sigEffect, n = 5))

cat("\nReturn last 5 rows of [merged_MA_sigEffect]\n")
print(tail(merged_MA_sigEffect, n = 5))

cat("\nSize of [merged_MA_sigEffect]: (# of rows, # of columns)\n")
print(dim(merged_MA_sigEffect))
cat("\n")


############################################################
### Identify SNPs with HetPVal < 5e-8
############################################################

merged_MA_sigHetPVal <- merged_MA %>%
  dplyr::filter(HetPVal < genome_sig)

cat("\nReturn first 5 rows of [merged_MA_sigHetPVal]\n")
print(head(merged_MA_sigHetPVal, n = 5))

cat("\nReturn last 5 rows of [merged_MA_sigHetPVal]\n")
print(tail(merged_MA_sigHetPVal, n = 5))

cat("\nSize of [merged_MA_sigHetPVal]: (# of rows, # of columns)\n")
print(dim(merged_MA_sigHetPVal))
cat("\n")


############################################################
### Get intersection of [significant effect] and [significant HetPVal SNPs]
############################################################

merged_MA_sigEffectHetPVal <- merged_MA_sigEffect %>%
  dplyr::filter(MarkerName %in% merged_MA_sigHetPVal$MarkerName)

cat("\nReturn first 5 rows of [merged_MA_sigEffectHetPVal]\n")
print(head(merged_MA_sigEffectHetPVal, n = 5))

cat("\nReturn last 5 rows of [merged_MA_sigEffectHetPVal]\n")
print(tail(merged_MA_sigEffectHetPVal, n = 5))

cat("\nSize of [merged_MA_sigEffectHetPVal]: (# of rows, # of columns)\n")
print(dim(merged_MA_sigEffectHetPVal))

cat("\n[chr]\n")
print(table(merged_MA_sigEffectHetPVal$chr))


############################################################
### Check how many significant SNPs fall into each FUMA Genomic Risk Locus
############################################################

genomicriskloci <- FUMA_results_GenomicRiskLoci %>% dplyr::filter(chr == CHR)

cat("\n[Identified Genomic Risk Loci]\n")
print(head(genomicriskloci, n = 5))

snps_to_check <- merged_MA_sigEffectHetPVal$position

genomicriskloci$SNP_count <- sapply(seq_len(nrow(genomicriskloci)), function(i) {
  sum(
    snps_to_check >= genomicriskloci$meta_start_hg19[i] &
      snps_to_check <= genomicriskloci$meta_end_hg19[i]
  )
})

cat("\n[Number of SNPs in merged_MA_sigEffectHetPVal in the identified Genomic Risk Loci]\n")
print(genomicriskloci)


############################################################
### Prepare data for effect-size flipping
############################################################

# Dynamic column names
col1 <- paste0("Effect_", group_1)
col2 <- paste0("Effect_", group_2)

# Use group 1 effect as reference direction
merged_MA$compare_effect <- merged_MA[[col1]]

# Flip signs so group 1 effect is always positive
merged_MA_flip <- merged_MA %>%
  mutate(
    !!sym(col1) := if_else(compare_effect < 0, -!!sym(col1), !!sym(col1)),
    !!sym(col2) := if_else(compare_effect < 0, -!!sym(col2), !!sym(col2))
  )

cat("\nReturn first 5 rows of [merged_MA]\n")
print(head(merged_MA, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA))

cat("\nReturn first 5 rows of [merged_MA_flip]\n")
print(head(merged_MA_flip, n = 5))
cat("\nSize of table: (# of rows, # of columns)\n")
print(dim(merged_MA_flip))


############################################################
### Keep only genomic risk loci with at least one SNP to highlight
############################################################

genomicriskloci_to_plot <- genomicriskloci %>% dplyr::filter(SNP_count > 0)

cat("\n[Genomic Risk Loci to Plot]\n")
print(genomicriskloci_to_plot)

if (nrow(genomicriskloci_to_plot) == 0) {
  cat("No genomic risk loci contain SNPs from merged_MA_sigEffectHetPVal on CHR", CHR, 
      ". Skipping this chromosome.\n")
  next
}


############################################################
### Loop over each Genomic Risk Locus and create separate plots
############################################################

for (i in seq_len(nrow(genomicriskloci_to_plot))) {

  locus_i <- genomicriskloci_to_plot[i, ]

  ############################################################
  ### Define locus ID
  ############################################################

    locus_id <- i
    locus_start_pos <- locus_i$meta_start_hg19
    locus_end_pos <- locus_i$meta_end_hg19

  cat("\n====================================================\n")
  cat("Plotting CHR", CHR, "Genomic Risk Locus", locus_id, "\n")
  cat("Locus start:", locus_start_pos, "\n")
  cat("Locus end:", locus_end_pos, "\n")
  cat("====================================================\n")

  ############################################################
  ### Get hg38 and hg 19 position of FUMA lead SNP for plot title
  ############################################################
  lead_snp <- as.character(locus_i$leadSNPs)
  lead_pos_hg19 <- locus_i$pos_hg19
  lead_info <- data_temp1 %>%
                    dplyr::filter(MarkerName == lead_snp) %>%
                    dplyr::slice(1)
  lead_chr <- lead_info$chr
  lead_pos_hg38 <- lead_info$position
    
  cat("\n====================================================\n")
  cat("Lead SNP:", lead_snp, "\n")
  cat("Lead SNP chromosome:", lead_chr, "\n")
  cat("Lead SNP position (hg38):", lead_pos_hg38, "\n")
  cat("Lead SNP position (hg19):", lead_pos_hg19, "\n")
  cat("====================================================\n")

  ############################################################
  ### Significant SNPs inside this specific locus
  ############################################################

  sig_snps_locus <- merged_MA_sigEffectHetPVal %>%
    dplyr::filter(
      position >= locus_i$meta_start_hg19,
      position <= locus_i$meta_end_hg19
    )

  snp_mark_locus <- sig_snps_locus$MarkerName
    
  cat("\nNumber of highlighted SNPs in this locus:", length(snp_mark_locus), "\n")

  if (length(snp_mark_locus) == 0) {
    cat("Skipping this locus because no SNPs are available to highlight.\n")
    next
  }
    

 ############################################################
 ### Define plotting window based on significant SNP positions
 ############################################################

 bp_min <- min(sig_snps_locus$position, na.rm = TRUE) - bp_extra
 bp_max <- max(sig_snps_locus$position, na.rm = TRUE) + bp_extra

 cat("\nPlotting window based on significant SNP positions:\n")
 cat("bp_min =", bp_min, "\n")
 cat("bp_max =", bp_max, "\n")
    
    
############################################################
### High-PIP SNPs inside this specific locus
############################################################

high_pip_snps_locus <- merged_MA_flip %>%
  dplyr::filter(
    MarkerName %in% data_high_pip_rsid$MarkerName,
    position >= bp_min,
    position <= bp_max
  )

high_pip_mark_locus <- high_pip_snps_locus$MarkerName

cat("\nNumber of high-PIP SNPs in this locus:",
    length(high_pip_mark_locus), "\n")
    

  ############################################################
  ### SNP with lowest HetPVal inside this specific locus
  ############################################################

  # SNP with the lowest HetPVal in this locus
  lowest_het_snp_locus <- sig_snps_locus %>%
    dplyr::slice_min(
      order_by = HetPVal,
      n = 1,
      with_ties = FALSE
    )%>%
    dplyr::mutate(
      GenomicRiskLocus = locus_id
    )
    
  cat("\nSNP with lowest HetPVal in this locus:\n")
  print(lowest_het_snp_locus)

    
  # add start and end position to [lowest_het_snp_locus]
  lowest_het_snp_locus_new <- lowest_het_snp_locus %>%
                                    mutate(meta_start_hg19 = locus_start_pos,
                                          meta_end_hg19 = locus_end_pos)

  # Save to full list
  lowest_het_snp_list[[paste0("chr", CHR, "_locus", locus_id)]] <- lowest_het_snp_locus_new

  ############################################################
  ### Filter SNPs to locus window
  ############################################################

  merged_MA_small <- merged_MA_flip %>%
    dplyr::filter(
      position >= bp_min,
      position <= bp_max
    )

  cat("\nReturn first 5 rows of [merged_MA_small]\n")
  print(head(merged_MA_small, n = 5))

  cat("\nReturn last 5 rows of [merged_MA_small]\n")
  print(tail(merged_MA_small, n = 5))

  cat("\nSize of [merged_MA_small]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small))
  cat("\n")


  ############################################################
  ### Filter SNPs to highlight
  ############################################################

  merged_MA_small_mark <- merged_MA_flip %>%
    dplyr::filter(MarkerName %in% snp_mark_locus)

  cat("\nReturn first 5 rows of [merged_MA_small_mark]\n")
  print(head(merged_MA_small_mark, n = 5))

  cat("\nReturn last 5 rows of [merged_MA_small_mark]\n")
  print(tail(merged_MA_small_mark, n = 5))

  cat("\nSize of [merged_MA_small_mark]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small_mark))
  cat("\n")


  ############################################################
  ### Create long-form data frame for Panel A
  ############################################################

  merged_MA_small_long <- merged_MA_small %>%
    dplyr::select(MarkerName, chr, position, !!sym(col1), !!sym(col2)) %>%
    pivot_longer(
      cols = c(!!sym(col1), !!sym(col2)),
      names_to = "group",
      values_to = "Effect"
    )

  cat("\nReturn first 5 rows of [merged_MA_small_long]\n")
  print(head(merged_MA_small_long, n = 5))

  cat("\nSize of [merged_MA_small_long]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small_long))


  ############################################################
  ### Mark SNPs to highlight
  ############################################################

  merged_MA_small_mark_long <- merged_MA_small_long %>%
    dplyr::filter(MarkerName %in% snp_mark_locus)

  cat("\nReturn first 5 rows of [merged_MA_small_mark_long]\n")
  print(head(merged_MA_small_mark_long, n = 5))

  cat("\nSize of [merged_MA_small_mark_long]: (# of rows, # of columns)\n")
  print(dim(merged_MA_small_mark_long))


  ############################################################
  ### Create high-PIP long data: Mark high-PIP SNPs
  ############################################################
  merged_MA_highPIP <- merged_MA_small %>%
    dplyr::filter(MarkerName %in% high_pip_mark_locus)

  merged_MA_highPIP_long <- merged_MA_small_long %>%
    dplyr::filter(MarkerName %in% high_pip_mark_locus)


  ############################################################
  ### Panel A: Marginal Genetic Effects Plot
  ############################################################

  effect_plot <- ggplot() +
    geom_segment(data = merged_MA_small,
                 aes(x = position,
                     xend = position,
                     y = !!sym(col1),
                     yend = !!sym(col2)),
                 linewidth = 0.25, alpha = 0.5, color = "gray75") +
    geom_point(data = dplyr::filter(merged_MA_small_long, group == col1),
               aes(position, Effect, color = paste0(Group_1, "_gray")),
               shape = 15, size = 2.3) +
    geom_point(data = dplyr::filter(merged_MA_small_long, group == col2),
               aes(position, Effect, color = paste0(Group_2, "_gray")),
               shape = 15, size = 2.3) +
    geom_segment(data = merged_MA_small_mark,
                 aes(x = position,
                     xend = position,
                     y = !!sym(col1),
                     yend = !!sym(col2)),
                 linewidth = 0.8,alpha = 0.5,color = "purple",show.legend = FALSE) +
    geom_hline(yintercept = 0,color = "gray",linewidth = 0.4,linetype = "dashed") +
    geom_point(data = dplyr::filter(merged_MA_small_mark_long, group %in% c(col1, col2)),
               aes(position, Effect, color = group),
               shape = 15,size = 2.2) +
    # High-PIP effect points
    geom_point(data = dplyr::filter(merged_MA_highPIP_long, group == col1),
               aes(position, Effect, color = "High-PIP: Group 1"),
               shape = 15, size = 2.8) +
    geom_point(data = dplyr::filter(merged_MA_highPIP_long, group == col2),
               aes(position, Effect, color = "High-PIP: Group 2"),
               shape = 15, size = 2.8) +    
    # Connect high-PIP SNP effects
    geom_segment(data = merged_MA_highPIP,
                 aes(x = position,
                     xend = position,
                     y = !!sym(col1),
                     yend = !!sym(col2)),
                 linewidth = 0.8, color = "#FFD600", show.legend = FALSE) +
    scale_color_manual(name = "",
                       values = c(
                       setNames("gray30", paste0(Group_1, "_gray")),
                       setNames("red", col1),
                       setNames("gray60", paste0(Group_2, "_gray")),
                       setNames("blue", col2),
                        "High-PIP: Group 1" = "#FFD600",
                        "High-PIP: Group 2" = "#FFD600"),
                      breaks = c(paste0(Group_1, "_gray"), col1,
                                 paste0(Group_2, "_gray"), col2,
                                 "High-PIP: Group 1"),
                      labels = c("",paste0(Group_1,"     "),
                                 "",paste0(Group_2,"     "),"High-PIP")) +
    guides(color = guide_legend(byrow = TRUE,nrow = 1,
                                keywidth = unit(0, "lines"))) +
    scale_x_continuous(labels = label_comma()) +
    coord_cartesian(xlim = c(bp_min, bp_max)) +
    labs(x = NULL, y = "Marginal Effects",
      title = paste0(Group_1, " vs ", Group_2, 
                     ", Chr", lead_chr, ":", bp_min, "-", bp_max)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 16),
      legend.spacing.x = unit(0, "cm"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.ticks.x = element_blank()
    )

  print(effect_plot)


  ############################################################
  ### Panel B: HetPVal Plot
  ############################################################

  df <- merged_MA_small %>%
    dplyr::mutate(
      HetPval = as.numeric(HetPVal),
      logHetP = -log10(HetPval)
    ) %>%
    dplyr::filter(!is.na(logHetP), is.finite(logHetP))

  hetpval_plot <- ggplot(df, aes(x = position, y = logHetP)) +
    geom_point(shape = 21,fill = "gray80",color = "black",stroke = 0.5,size = 2) +
    geom_hline(yintercept = -log10(genome_sig),
               linetype = "dashed",color = "gray60",linewidth = 0.4) +
    scale_x_continuous(labels = label_comma()) +
    coord_cartesian(xlim = c(bp_min, bp_max)) +
    labs(x = NULL,
         y = expression(-log[10]("HetP-value")),
         title = "") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.ticks.x = element_blank()
    ) +
    geom_point(data = df %>% dplyr::filter(MarkerName %in% snp_mark_locus),
      aes(position, logHetP),
      shape = 21,fill = "deepskyblue",color = "black",stroke = 0.5,size = 2
    )

  print(hetpval_plot)


  ############################################################
  ### Panel C: Mapped Genes
  ############################################################

  data_locus <- FUMA_results %>%
    dplyr::filter(chr == CHR) %>%
    dplyr::filter(type == "protein_coding") %>%
    dplyr::filter(
      end_gene >= bp_min,
      start_gene <= bp_max
    ) %>%
    mutate(
      start_mb = start_gene,
      end_mb = end_gene,
      plot_start = pmax(start_mb, bp_min),
      plot_end   = pmin(end_mb, bp_max),
      label_x = (plot_start + plot_end) / 2
    ) %>%
    arrange(start_mb) %>%
    mutate(y = row_number())

  data_name <- paste0("genes_chr", CHR, "_locus", locus_id)
  assign(data_name, data_locus)

  cat("\nReturn first 5 rows of mapped genes for this locus\n")
  print(head(get(data_name), n = 5))

  cat("\nReturn last 5 rows of mapped genes for this locus\n")
  print(tail(get(data_name), n = 5))

  cat("\nSize of mapped genes table: (# of rows, # of columns)\n")
  print(dim(get(data_name)))


  if (nrow(data_locus) > 0) {

    gene_track_chr <- ggplot(data_locus) +
      geom_segment(aes(x = start_mb,xend = end_mb,y = y,yend = y),
                   linewidth = 1.2,color = "black") +
      geom_text(aes(x = label_x,y = y,label = symbol),
                    size = 4,vjust = -0.6) +
      scale_y_continuous(breaks = NULL) +
      scale_x_continuous(labels = label_comma()) +
      coord_cartesian(
        xlim = c(bp_min, bp_max),
        ylim = c(min(data_locus$y) - 0.5, max(data_locus$y) + 0.5)) +
      labs(x = paste0("Base Pair Position on Chromosome ", CHR),
           y = "Mapped Genes") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 8)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)
      )

  } else {


    gene_track_chr <- ggplot() +
      coord_cartesian(xlim = c(bp_min, bp_max)) +
      scale_x_continuous(labels = label_comma()) +
      labs(x = paste0("Base Pair Position on Chromosome ", CHR),
           y = "Mapped Genes") +
      annotate("text",
                x = mean(c(bp_min, bp_max)),
                y = 1,
                label = "No protein-coding genes in this window",size = 5) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 8)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  print(gene_track_chr)


  ############################################################
  ### Combine Panels A, B, and C
  ############################################################

  combined_plot <- effect_plot / hetpval_plot / gene_track_chr +
    plot_layout(heights = c(1, 1, 1))

  print(combined_plot)


  ############################################################
  ### Save combined plot
  ############################################################

  output_file <- paste0(root_dir, "/", type,"/Output/Figures/locuszoom_chr", CHR,
                        "_locus", locus_id,"_", type2,"_sigHetPVal_FUMA_zoomedin.png")

  ggsave(output_file,combined_plot,width = 10,height = 9,dpi = 500)

  cat("\nSaved plot:\n")
  cat(output_file, "\n")
}
    
cat("\n*-----------------------------------------------------------------------------------*\n")

} # for: CHR


############################################################
### Combine lowest-HetPVal SNPs across all loci
############################################################

lowest_het_snp_all <- dplyr::bind_rows(
  lowest_het_snp_list,
  .id = "locus_list"
)

lowest_het_snp_all <- lowest_het_snp_all %>%
 dplyr::mutate(type=type2) %>% 
 dplyr::arrange(chr, position)

cat("\n[Lowest HetPVal SNP from each Genomic Risk Locus]\n")
print(lowest_het_snp_all)


# Save as txt
file_path = paste0(root_dir, "/", type,"/Output/", type3, "_Lowest_HetPval")
write.table(lowest_het_snp_all,
      file = paste0(file_path, "/LowestHetPValSNP_",type2,".txt"), 
      sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE
     )
