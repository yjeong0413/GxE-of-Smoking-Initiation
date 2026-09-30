R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#ANCESTRY="european"

cat("Set directories\n")
start_dir="/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024"
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
 #type2: defined below
type3 <- "Generation"
category <- "generation"
 #subgroups <- paste0("generation",1:4)
subgroups <- paste0("birthyear_bin",4:19)
subgroups

for (SUBGROUP in subgroups) {
    
    cat("\n")
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))

    file_path <- paste0(root_dir,"/",
                        type,"/METAL_cross-ancestry_",SUBGROUP,"_allancestry_MAresults_cleaned2.txt")
                               
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

# Read MA results
file_path <- paste0(root_dir, "/",type,"/meta_locus_table_for_locuszoom_with_original_hg19.tsv")
cat(paste0("\n1. [", file_path, "]\n"))
cat("1) Read the data.\n")
   data_temp <- fread(file_path, sep = "\t", header = T)

   data <- data_temp %>%
                        filter(locus_class == paste0(category,"_main"))

    # assign data to [data_meta_loci]
    data_name <- paste0("data_meta_loci")
    assign(data_name,data, envir = .GlobalEnv)

    cat("\n[",data_name,"]\n")
    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

   # pick necessary columns
   data_simple <- data %>%
                    select(chr_hg38, meta_start_hg38, meta_end_hg38) %>%
                    filter(!is.na(chr_hg38),
                          !is.na(meta_start_hg38),
                          !is.na(meta_end_hg38),) %>%
                    distinct()

    # assign data to [data_meta_loci_simple]
    data_name <- paste0("data_meta_loci_simple")
    assign(data_name,data_simple, envir = .GlobalEnv)

    cat("\n[",data_name,"]\n")
    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 10))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(get(data_name), n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

cat("\nNumber of unique main-effect loci identified using [chr_hg38, meta_start_hg38, meta_end_hg38]:\n")
 print(get(data_name) %>%
       distinct(chr_hg38, meta_start_hg38, meta_end_hg38) %>%
       nrow()
      )


file_path <- paste0(root_dir,"/",
                    type,"/METAL_cross-ancestryMA_",type,"_MAresults_cleaned.txt")
                               
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
    data_name <- paste0("cleaned_birthyearMA")
    assign(paste0("cleaned_birthyearMA"), cleaned_data, envir = .GlobalEnv)
        
    cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))


############################################################
### Find SNP with lowest HetPVal within each meta locus
############################################################

library(data.table)

# Convert to data.table
setDT(data_meta_loci_simple)
setDT(cleaned_birthyearMA)

############################################################
### Add unique ID for each meta locus
############################################################
data_meta_loci_simple[, locus_id := .I]

############################################################
### Match SNPs falling within each meta locus
############################################################
SNPs_in_meta_loci <- cleaned_birthyearMA[
  data_meta_loci_simple,
  on = .(
    chr = chr_hg38,
    position >= meta_start_hg38,
    position <= meta_end_hg38
  ),
  nomatch = 0,
  allow.cartesian = TRUE,
  .(
    locus_id = i.locus_id,
    chr_hg38 = i.chr_hg38,
    meta_start_hg38 = i.meta_start_hg38,
    meta_end_hg38 = i.meta_end_hg38,
    MarkerName,
    position = x.position,
    `P-value`,
    HetPVal,
    Effect,
    StdErr
  )
]

############################################################
### Select SNP with lowest PVal within each meta locus
############################################################
SNP_lowest_PVal <- SNPs_in_meta_loci[
  !is.na(`P-value`)
][
  order(locus_id, `P-value`)
][
  , .SD[1], by = locus_id
]

############################################################
### Check results
############################################################
cat("\nNumber of meta loci:\n")
 print(nrow(data_meta_loci_simple))
cat("\nNumber of representative SNPs:\n")
 print(nrow(SNP_lowest_PVal))
cat("\nFirst 10 representative SNPs:\n")
 print(head(SNP_lowest_PVal, 10))

############################################################
### Assign data
############################################################
SNP_mainmod_list <- copy(SNP_lowest_PVal)
cat("\nNumber of representative SNPs:\n")
 print(nrow(SNP_lowest_PVal))


############################################################
############################################################
### Birthyear bins 4-19: Main-effect 1,753 SNPs
############################################################
############################################################

library(data.table)
library(dplyr)
library(nlme)
library(ggplot2)

############################################################
### 0. Settings
############################################################
birthyear_bins <- 4:19

# Name used for output files
type2 <- "birthyear_bin4_19"

# Output directory
file_path <- paste0(root_dir, "/", type, "/Output/MainEffect_1753_Birthyear_Bin4_19")

dir.create(file_path, recursive = TRUE, showWarnings = FALSE)

############################################################
############################################################
### 1. Main-effect SNP list
############################################################
############################################################
main_effect_snps <- unique(SNP_mainmod_list$MarkerName)

cat("\nNumber of unique main-effect SNPs:\n")
print(length(main_effect_snps))

############################################################
############################################################
### 2. Combine birthyear_bin4 - birthyear_bin19
###    and keep only the 1,753 main-effect SNPs
############################################################
############################################################
raw_effect_birthyear <- bind_rows(
  lapply(birthyear_bins, function(i) {
    data_name <- paste0("cleaned_birthyear_bin", i)
    cat("\nReading object: ", data_name, "\n")
    dat <- get(data_name)
    dat %>%
      filter(MarkerName %in% main_effect_snps) %>%
      mutate(
        birthyear_bin = i,
        birthyear_num = i - 3
      ) #mutate
  }) #lapply

) #bind_rows

############################################################
### 3. Harmonize risk-increasing direction
###
### If average effect across the 16 birthyear bins < 0,
### flip the effect sign and allele frequency.
############################################################

raw_effect_birthyear <- raw_effect_birthyear %>%
                            group_by(MarkerName) %>%
                                mutate(
                                    # Average effect across birthyear bins 4-19
                                    mean_effect = mean(Effect_new, na.rm = TRUE),
                                    # Whether this SNP should be flipped
                                    flip_effect = mean_effect < 0,
                                    # Harmonized effect size
                                    Effect_harmonized = if_else(flip_effect,
                                                                -Effect_new, Effect_new),
                                    # Harmonized effect allele frequency
                                    Freq1_harmonized = if_else(flip_effect,
                                                               1 - Freq1_new, Freq1_new)) %>% 
                            ungroup()

############################################################
### 4. Check harmonization
############################################################
harmonization_check <- raw_effect_birthyear %>%
                            distinct(MarkerName, mean_effect, flip_effect)

cat("\nNumber of SNPs flipped vs. not flipped:\n")
print(table(harmonization_check$flip_effect))

############################################################
### 5. Basic checks
############################################################
cat("\n---------------------------------------------\n")
cat("Number of observations by birthyear bin\n")
cat("---------------------------------------------\n")
print(table(raw_effect_birthyear$birthyear_bin))

cat("\nNumber of unique SNPs:\n")
print(length(unique(raw_effect_birthyear$MarkerName)))

cat("\nDimension of combined dataset:\n")
print(dim(raw_effect_birthyear))

# If every one of the 1,753 SNPs occurs in all 16 bins:
# 1753 * 16 = 28,048 rows

############################################################
### 6. Check number of bins available for each SNP
############################################################
snp_check <- raw_effect_birthyear %>%
  group_by(MarkerName) %>%
  summarise(
    n_birthyear_bins = n_distinct(birthyear_bin),
    .groups = "drop"
  ) #summarise

cat("\nNumber of birthyear bins available per SNP:\n")
print(table(snp_check$n_birthyear_bins))

############################################################
### 7. Keep SNPs observed in all 16 bins
############################################################
complete_snps <- snp_check %>%
  filter(n_birthyear_bins == 16) %>%
  pull(MarkerName)

cat("\nNumber of SNPs observed in all 16 bins:\n")
print(length(complete_snps))
cat("\nFirst 5 rows:\n")
print(head(complete_snps, n=5))

raw_effect_birthyear <- raw_effect_birthyear %>% 
                            filter(MarkerName %in% complete_snps)
cat("\nDimension:\n")
print(dim(raw_effect_birthyear))
cat("\nFirst 5 rows:\n")
print(head(raw_effect_birthyear, n=5))

############################################################
### 8. Compute 95% CI for effect sizes
############################################################
raw_effect_birthyear <- raw_effect_birthyear %>%
    mutate(
        CI_lower_harmonized = Effect_harmonized - 1.96 * StdErr,
        CI_upper_harmonized = Effect_harmonized + 1.96 * StdErr
    ) #mutate

############################################################
### 9. Sort data
############################################################
raw_effect_birthyear <- raw_effect_birthyear %>%
    arrange(MarkerName, birthyear_num)

cat("\nDimension:\n")
print(dim(raw_effect_birthyear))
cat("\nFirst 5 rows:\n")
print(head(raw_effect_birthyear, n=5))

############################################################
### 10. Final checks
############################################################
cat("\nFinal number of SNPs:\n")
print(length(unique(raw_effect_birthyear$MarkerName)))

cat("\nFinal number of observations:\n")
print(nrow(raw_effect_birthyear))

# If every one of the 1,670 SNPs occurs in all 16 bins:
# 1670 * 16 = 26,720 rows

############################################################
############################################################
### 11. Initialize results table
############################################################
############################################################
results <- data.frame(
  snp = character(),
  chr = numeric(),
  position = numeric(),
  eff_intercept = numeric(),
  eff_slope = numeric(),
  eff_pvalue = numeric(),
  eaf_intercept = numeric(),
  eaf_slope = numeric(),
  eaf_pvalue = numeric(),
  stringsAsFactors = FALSE
) #data.frame

############################################################
############################################################
### 12. Split data by SNP
############################################################
############################################################
plots <- split(
  raw_effect_birthyear,
  raw_effect_birthyear$MarkerName
)

############################################################
############################################################
### 13. Run GLS separately for each SNP
############################################################
############################################################
cnt = 0

for (snp in names(plots)) {

  cnt = cnt+1
  cat("\n[Plot number]:\n")
  print(cnt)
  ##########################################################
  ### Data for one SNP
  ##########################################################

  d <- plots[[snp]] %>% arrange(birthyear_num)

   ##########################################################
   ### SNP information 
   ##########################################################
   chr <- d$chr[1]
   position <- d$position[1]

   cat("\n")
   cat("====================================================\n")
   cat(paste0("SNP: ", snp, "\n"))
   cat(paste0("Chr: ", chr, "\n"))
   cat(paste0("Position: ", position, "\n"))
   cat("====================================================\n")

##########################################################
  ##########################################################
  ### A. EFFECT SIZE
  ##########################################################
  ##########################################################

  ##########################################################
  ### A1. GLS model
  ##########################################################

  model_eff <- gls(
    Effect_harmonized ~ birthyear_num,
    data = d,
    weights = varFixed(~ StdErr^2)
  ) #gls

  cat("\n")
  cat(paste0("[",snp,": Effect Size]\n"))

  print(summary(model_eff))

  ##########################################################
  ### A2. Extract coefficients
  ##########################################################

  eff_intercept <- coef(model_eff)[1]

  eff_slope <- coef(model_eff)[2]

  eff_pvalue <- summary(model_eff)$tTable["birthyear_num","p-value"]

  ##########################################################
  ### A3. P-value label
  ##########################################################

  eff_p_label <- ifelse(
      eff_pvalue < 0.001,
      "p-value for trend < 0.001",
    paste0("p-value for trend = ",formatC(eff_pvalue, format = "f", digits = 3))
  ) #ifelse

  ##########################################################
  ### A4. Fitted values
  ##########################################################

  d$fitted_eff <- predict(model_eff, newdata = d)

  ##########################################################
  ### A5. Effect-size plot
  ##########################################################

  p_eff <- ggplot(d, aes(x = birthyear_num)) +
    ########################################################
    ### 95% CI
    ########################################################
    geom_ribbon(aes(ymin = CI_lower_harmonized, ymax = CI_upper_harmonized, fill = "95% CI"), 
                alpha = 0.25) +
    ########################################################
    ### Raw effect trajectory
    ########################################################
    geom_line(aes(y = Effect_harmonized, color = "Raw"), linewidth = 1) +
    geom_point(aes(y = Effect_harmonized, color = "Raw"), size = 3) +
    ########################################################
    ### Fitted GLS trend
    ########################################################
    geom_line(aes(y = fitted_eff, color = "Fitted"), linewidth = 1) +
    ########################################################
    ### Zero reference line
    ########################################################
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ########################################################
    ### X-axis
    ########################################################
    scale_x_continuous(breaks = 1:16, labels = paste0("bin", 4:19)) +
    ########################################################
    ### Labels
    ########################################################
    labs(title = paste0("Effect Size (", snp, ", Birthyear Bins 4-19)"), 
        x = NULL, y = "Effect Size", color = "", fill = "") +
    scale_color_manual(values = c("Raw" = "black", "Fitted" = "blue")) +
    scale_fill_manual(values = c("95% CI" = "gray60")) +
    ########################################################
    ### P-value
    ########################################################
    annotate("text", x = Inf, y = Inf,
        label = eff_p_label, hjust = 1.15, vjust = 1.5, size = 5) +
    ########################################################
    ### Theme
    ########################################################
    theme_minimal(base_size = 18) +
    theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "top",
          legend.text = element_text(size = 13)
         ) #theme

  print(p_eff)
    
  ##########################################################
  ### A6. Save effect-size plot
  ##########################################################

  ggsave(filename = paste0(file_path, "/", type2, "_", snp, "_Raw_GLS_Effect_Size_main.png"), 
         plot = p_eff, width = 9, height = 6, dpi = 500)


  ##########################################################
  ##########################################################
  ### B. EFFECT ALLELE FREQUENCY
  ##########################################################
  ##########################################################
    
  ##########################################################
  ### B1. Compute SE and 95% CI for allele frequency
  ##########################################################
  d$Freq1_SE_harmonized <- sqrt(d$Freq1_harmonized * (1 - d$Freq1_harmonized)/(2 * d$N_tot))

  d$Freq1_CI_lower_harmonized <- d$Freq1_harmonized - 1.96 * d$Freq1_SE_harmonized
  d$Freq1_CI_upper_harmonized <- d$Freq1_harmonized + 1.96 * d$Freq1_SE_harmonized

  ##########################################################
  ### Keep frequency CI within 0-1
  ##########################################################
  d$Freq1_CI_lower_harmonized <- pmax(d$Freq1_CI_lower_harmonized, 0)
  d$Freq1_CI_upper_harmonized <- pmin(d$Freq1_CI_upper_harmonized, 1)

  ##########################################################
  ### B2. GLS model for allele frequency
  ##########################################################
  model_freq <- gls(
    Freq1_harmonized ~ birthyear_num,
    data = d,
    weights = varFixed(~ Freq1_harmonized * (1 - Freq1_harmonized)/(2 * N_tot))
  ) #gls

  cat("\n")
  cat(paste0("[", snp, ": Effect Allele Frequency]\n"))
  print(summary(model_freq))

  ##########################################################
  ### B3. Extract coefficients
  ##########################################################
  eaf_intercept <- coef(model_freq)[1]
  eaf_slope <- coef(model_freq)[2]
  eaf_pvalue <- summary(model_freq)$tTable["birthyear_num", "p-value"]

  ##########################################################
  ### B4. P-value label
  ##########################################################
  eaf_p_label <- ifelse(eaf_pvalue < 0.001,
                        "p-value for trend < 0.001",
                        paste0("p-value for trend = ", formatC(eaf_pvalue, format = "f", digits = 3))
                       ) #ifelse

  ##########################################################
  ### B5. Fitted values
  ##########################################################
  d$fitted_freq <- predict(model_freq, newdata = d)

  ##########################################################
  ### B6. Frequency plot
  ##########################################################
  p_freq <- ggplot(d, aes(x = birthyear_num)) +
    ########################################################
    ### 95% CI
    ########################################################
    geom_ribbon(aes(ymin = Freq1_CI_lower_harmonized, ymax = Freq1_CI_upper_harmonized, fill = "95% CI"), 
                alpha = 0.25) +
    ########################################################
    ### Raw EAF trajectory
    ########################################################
    geom_line(aes(y = Freq1_harmonized, color = "Raw"), linewidth = 1) +
    geom_point(aes(y = Freq1_harmonized, color = "Raw"), size = 3) +
    ########################################################
    ### Fitted GLS trend
    ########################################################
    geom_line(aes(y = fitted_freq, color = "Fitted"), linewidth = 1) +
    ########################################################
    ### X-axis
    ########################################################
    scale_x_continuous(breaks = 1:16, labels = paste0("bin", 4:19)) +
    ########################################################
    ### Labels
    ########################################################
    labs(title = paste0("Effect Allele Frequency (", snp, ", Birthyear Bins 4-19)"), 
         x = NULL, y = "Effect Allele Frequency", color = "", fill = "") +
    scale_color_manual(values = c("Raw" = "black", "Fitted" = "blue")) +
    scale_fill_manual(values = c("95% CI" = "gray60")) +
    ########################################################
    ### P-value
    ########################################################
    annotate("text", x = Inf, y = Inf, label = eaf_p_label, hjust = 1.15, vjust = 1.5, size = 5) +
    ########################################################
    ### Theme
    ########################################################
    theme_minimal(base_size = 18) +
    theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "top",
          legend.text = element_text(size = 13)
         ) #theme

  print(p_freq)

  ##########################################################
  ### B7. Save frequency plot
  ##########################################################
  ggsave(filename = paste0(file_path, "/", type2, "_", snp, "_Raw_GLS_Effect_Allele_Frequency_main.png"), 
         plot = p_freq, width = 9, height = 6, dpi = 500)

  ##########################################################
  ##########################################################
  ### C. Store results
  ##########################################################
  ##########################################################
  results <- rbind(results,
                   data.frame(
                       snp = snp,
                       chr = chr,
                       position = position,
                       eff_intercept = as.numeric(eff_intercept),
                       eff_slope = as.numeric(eff_slope),
                       eff_pvalue = as.numeric(eff_pvalue),
                       eaf_intercept = as.numeric(eaf_intercept),
                       eaf_slope = as.numeric(eaf_slope),
                       eaf_pvalue = as.numeric(eaf_pvalue))
                  ) #rbind

} # for: snp


############################################################
############################################################
### 14. Add genomic-locus information
############################################################
############################################################
SNP_mainmod_unique <- SNP_mainmod_list %>%
                            distinct(MarkerName, .keep_all=TRUE)
results <- results %>%
            left_join(
                SNP_mainmod_unique %>%
                select(MarkerName) %>% distinct(), 
                by = c("snp" = "MarkerName")
            ) # left_join

############################################################
############################################################
### 15. Save GLS results
############################################################
############################################################
results_file <- paste0(file_path, "/MainEffect_1753_Birthyear_Bin4_19_GLS_results.tsv")
fwrite( results, results_file, sep = "\t", quote = FALSE, na = "NA")

############################################################
### 16. Print final results
############################################################
cat("\n\n============================================\n")
cat("GLS analysis completed\n")
cat("============================================\n")

cat("\nNumber of SNPs analyzed:\n")
num_results <- nrow(results)
print(num_results)

cat("\nFirst 10 results:\n")
print(head(results, 10))

cat("\nResults saved to:\n")
cat(results_file, "\n")

############################################################
### 17. Print some statistics
############################################################
cat("\n\n============================================\n")
cat("Some statistics:\n")
cat("============================================\n")
 ############################################################
 ### Max, min, and mean of effect-size trends (slope)
 ############################################################
 cat("\n Max of effect-size trends (slope):\n")
 print(max(results$eff_slope))
 cat("\n Min of effect-size trends (slope):\n")
 print(min(results$eff_slope))
 cat("\n Mean of effect-size trends (slope):\n")
 print(mean(results$eff_slope))
 ############################################################
 ### Max, min, and mean of effect-size trends p-value
 ############################################################
 cat("\n Max of effect-size trends p-value:\n")
 print(max(results$eff_pvalue))
 cat("\n Min of effect-size trends p-value:\n")
 print(min(results$eff_pvalue))
 cat("\n Mean of effect-size trends p-value:\n")
 print(mean(results$eff_pvalue))
 ############################################################
 ### Number of significant effect-size trends (p < 0.05)
 ############################################################
 threshold_p <- 0.05
 cat("\nTotal number of SNPs tested:\n")
 print(num_results)
 cat("\nNumber of significant effect-size trends (p < ",threshold_p,"):\n")
 sum(results$eff_pvalue < threshold_p, na.rm = TRUE)
 cat("\nPercent significant:\n")
 print(sum(results$eff_pvalue < threshold_p, na.rm = TRUE)/num_results)

 cat("\nNumber of significant effect-size trends (p < ",threshold_p,") by direction:\n")
 results %>% filter(eff_pvalue < threshold_p) %>% summarise(total_significant = n(),
                                                            positive_slope = sum(eff_slope > 0),
                                                            negative_slope = sum(eff_slope < 0))

 cat("\nNumber of significant effect-size trends (p < ",threshold_p/num_results,") after multiple testing correction:\n")
 sum(results$eff_pvalue < threshold_p/num_results, na.rm = TRUE)


threshold_p <- 0.05
results %>%
  summarise(
    total = n(),
    mean_slope = mean(eff_slope, na.rm = TRUE),
    median_slope = median(eff_slope, na.rm = TRUE),

    slope_positive = sum(eff_slope > 0, na.rm = TRUE),
    slope_negative = sum(eff_slope < 0, na.rm = TRUE),

    slope_significant = sum(eff_pvalue < threshold_p, na.rm = TRUE),
    slope_sig_positive = sum(eff_pvalue < threshold_p & eff_slope > 0, na.rm = TRUE),
    slope_sig_negative = sum(eff_pvalue < threshold_p & eff_slope < 0, na.rm = TRUE),

    slope_sig_bonf = sum(eff_pvalue < threshold_p/total, na.rm = TRUE),
    slope_sig_pos_bonf = sum(eff_pvalue < threshold_p/total & eff_slope > 0, na.rm = TRUE),
    slope_sig_neg_bonf = sum(eff_pvalue < threshold_p/total & eff_slope < 0, na.rm = TRUE)
  ) #summarize

############################################################
### Histogram of GLS effect-size slopes
############################################################
# Note: The histogram itself does not tell us whether individual trends are statistically significant.
p_slope_hist <- ggplot(results, aes(x = eff_slope)) +
  geom_histogram(bins = 40, color = "black", fill = "gray70") +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Distribution of Effect-Size Trends Across Main-Effect SNPs",
       x = "GLS Slope of Effect Size",
       y = "Number of SNPs") +
  theme_minimal(base_size = 18) +
  theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
       ) #theme

print(p_slope_hist)

############################################################
### Save plot
############################################################
ggsave(filename = paste0(file_path, "/MainEffect_1753_Effect_Size_GLS_Slope_Histogram.png"), 
       plot = p_slope_hist, width = 8, height = 6, dpi = 500)

############################################################
### Slope distribution with significant SNPs overlaid
############################################################
#- Keep the histogram or density plot.
#- Overlay the slope of SNPs with [p < 0.05] as little ticks along the bottom.
#- This immeidately shows whether significant trends are concentrated in the negative tail, positive tail, or both
p_slope_sig_hist <- ggplot(results, aes(x = eff_slope)) +
    geom_histogram(bins = 40, color = "black", fill = "gray70") +
    geom_rug(data = results %>% filter(eff_pvalue < 0.05), sides = "b", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "Effect-Size Trend (GLS Slope)", y = "Density") +
    theme_minimal(base_size = 16) +
    theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
         ) #theme

print(p_slope_sig_hist)

############################################################
### Save plot
############################################################
ggsave(filename = paste0(file_path, "/MainEffect_1753_Effect_Size_GLS_Slope_Sig_Histogram.png"), 
       plot = p_slope_sig_hist, width = 8, height = 6, dpi = 500)

############################################################
### Volcano-style slope plot
############################################################
#- Put [eff_slope] on the x-axis and [-log10(eff_pvalue)] on the y-axis.
#- SNPs on the left have decreasing effects across birth cohorts.
#- SNPs on the right have increasing effects across birth cohorts.
#- SNPs toward the top have stronger evidence for a temporal trend.
#- Add a horizontal line at [p=0.05].
#- This combines [direction] + [magnitude] + [statistical significance] in one figure

# Data
results <- results %>% 
            mutate(neglog10_p = -log10(eff_pvalue),
                   significant = eff_pvalue < 0.05,
                   significant_bon = eff_pvalue < 0.05/num_results
                  ) #mutate

# Nominal threshold
threshold_p <- 0.05
nom_neglog10 <- -log10(threshold_p)
nom_neglog10
# Bonferroni threshold
bon_neglog10 <- -log10(threshold_p / num_results)
bon_neglog10

# Plot
p_slope_sig_scatter <- ggplot(results, aes(x = eff_slope, y = neglog10_p)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
     scale_color_manual(name = "Significance",
                        values = c("FALSE" = "gray10", "TRUE" = "red"),
                        labels = c("Not Significant", "Significant (p < 0.05)")) +
  #Slope = 0
  geom_vline(aes(xintercept = 0, linetype = "Slope = 0"), color = "black", 
             linewidth = 0.4) +
  #Nominal significance threshold (p=0.05)
  geom_hline(aes(yintercept = nom_neglog10, linetype = "-log10(p) = 1.301"), color = "blue",
             linewidth = 0.4) +
  #Bonferroni significance threhold
  geom_hline(aes(yintercept = bon_neglog10, linetype = "Bonferroni"), color = "green", 
             linewidth = 0.4) +
  scale_linetype_manual(name = "Reference Lines",
                        values = c("Slope = 0" = "solid", 
                                   "-log10(p) = 1.301" = "solid",
                                   "Bonferroni" = "solid"),
                        breaks = c("Slope = 0","-log10(p) = 1.301","Bonferroni"),
                        labels = c("Slope = 0",
                                   paste0("Nominal: -log10(p) = ", round(nom_neglog10, 4)),
                                   paste0("Bonferroni: -log10(p) = ", round(bon_neglog10, 4))),
                        guide = guide_legend(override.aes = list(color = c("black","blue","green"),
                                                                 linewidth = 0.4))) +
  labs(x = "Effect-Size Trend (GLS Slope)",
       y = expression(-log[10](p))) +
  theme_minimal(base_size = 16) +
  theme(panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
       ) #theme

print(p_slope_sig_scatter)

############################################################
### Save plot
############################################################
ggsave(filename = paste0(file_path, "/MainEffect_1753_Effect_Size_GLS_Slope_Sig_Scatter.png"), 
       plot = p_slope_sig_scatter, width = 8, height = 6, dpi = 500)








