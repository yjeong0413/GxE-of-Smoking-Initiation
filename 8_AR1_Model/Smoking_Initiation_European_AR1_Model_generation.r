R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo() 

# Define ANCESTRY & SUBGROUP
ANCESTRY="european"

#Define directories 
cat("Set directories\n")
org_dir="/YOUR PATH HERE/"

 # Smoking Initiation
 root_dir=paste0(org_dir,"/",ANCESTRY,"/CLEAN/GWAS_QC")
 output_dir=paste0(root_dir,"/AR1_Model/results")

cat("\n-------------------------------------------------------------------------------\n")

cat("\nLoad Libraries\n")
cat("\n-----------------------------------------------------------------------------------\n")
cat("*Load [data.table]\n")
library(data.table)
cat("\n")
cat("*Load [readxl]\n")
library(readxl)
cat("\n")
cat("*Load [dplyr]\n")
library(dplyr)
cat("\n")
cat("*Load [tidyverse]\n")
library(tidyverse)
cat("\n")
cat("*Load [purrr]\n")
library(purrr)
cat("\n")
cat("*Load [qqman]\n")
library(qqman)
cat("\n")
cat("\n*Load [ggplot2]\n")
library(ggplot2)
cat("\n")
cat("\n*Load [patchwork]\n")
library(patchwork)
cat("\n")
cat("\n*Load [gridExtra]\n")
library(gridExtra)
cat("\n")
cat("\n*Load [scales]\n")
library(scales)
cat("\n")
cat("\n*Load [xtable]\n")
library(xtable)
cat("\n")
cat("\n*Load [TwoSampleMR]\n")
#install.packages("TwoSampleMR")
library(TwoSampleMR)
cat("\n*Load [ieugwasr]\n")
library(ieugwasr)
 # Regenerate token here [https://api.opengwas.io/] if clumping does not work.
 Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJqZW9uZzIxMkBwdXJkdWUuZWR1IiwiaWF0IjoxNzcxODcyNzAwLCJleHAiOjE3NzMwODIzMDB9.O03JhAO0BVNYz-xgkbTBIvm3BcEj8NxmDe_EOIVsma_j8hX_26EdulXd1btiXXGCYOMZKm_PhiQFarSSIDxoq-tOhXA_8GJZgtfpxZovL3LBl-e2Srkdk4y49nNuPYUKigIyg07BeTVJALsteS1kr2yheY6IKTyiTD2eN4CKlwRhZdLUWcFT2ssQeILDm6dkunLF5mObtcOac_57g9i9OXQiWNReGdUYJq9LmUmQqPan0IBfR9G7FSXelIr3p3TVziAsQ0B-ehDA--8HIqf8Fusw6QiNKR17156u34LBUMfvurTwn2ME4t47h9tGnD6lLvCbMKF6ngTssbDlPw2Pgw")
cat("\n")
cat("\n*Load [devtools]\n")
library(devtools)
cat("\n")
cat("\n*Load [GenomicSEM]\n")
#install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)
cat("\n")
cat("\n*Load [corrplot]\n")
#install.packages("corrplot")
library(corrplot)
cat("\n-----------------------------------------------------------------------------------\n")


# Define subgroup
type <- "birthyear"
type_short <- "gen"
category <- "groupingBIN"
subgroups <- paste0(category,1:4)
subgroups

# --------------------------------
# Dictionary (GWAS -> phenotypes)
# --------------------------------

# (03/24/2026) To above V > 1000 warning, I think it will be better if we group phenotypes by 
# sociodemo, health conditions, substance use. Then, run the model.

#define [GWAS_dict] with [Names=GWAS] and [Values=phenotypes]
GWAS_dict <- list(
  Okbay_et_al_2022  = c("EA"),  
  Becker_et_al_2021 = c("BMI", "HEIGHT", "CP", "AFB", "AUDIT", "ASTHMA",
                        "ASTECZRHI", "CANNABIS", "COPD", "DEPSYM", "ECZEMA",
                        "HAYFEVER", "MIGRAINE", "NEARSIGHTED", "SELFHEALTH",
                        "FAMSAT", "FINSAT", "FRIENDSAT", "WORKSAT", "LONELY",
                        "NEURO", "RELIGATT", "RISK", "SWB"),
  Gupta_et_al_2024  = c("NEURO", "AGREE", "CONSC", "OPEN", "EXTRA"),
  Liu_et_al_2019    = c("SI", "SC", "AGESI", "CPD", "DPW"),
  NEALE_v2          = c("TDI", "BMI", "ASTHMA", "COPD", "HAYFEVER", "HDL", 
                        "DIABII", "STROKE", "HAMI", "LUNG", "HYPTENS", "LONELY", 
                        "DEP", "RELIGATT"),
  PGC               = c("SCHIZ", "ADHD", "MDD", "ASD", "BIP")

) #list: GWAS_dict (Names=GWAS, Values=phenotypes)

# --------------------------------
# Generate separate trait dictionary for [Substance Use], [Sociodemographics], [Psychiatric and Medical Conditions] 
# --------------------------------
# define subset function
subset_GWAS_pairs <- function(GWAS_dict, pair_list) {
    
    result <- list()
    
    for (gwas_name in names(pair_list)) {
        if (gwas_name %in% names(GWAS_dict)) {
            result[[gwas_name]] <- intersect(GWAS_dict[[gwas_name]], pair_list[[gwas_name]])
        } #if
    } #for
    
    # remove empty entries
    result[sapply(result, length) > 0]
    
} #function: subset_GWAS_pairs

# -------------------------
# [Substance Use]
# -------------------------
substance_pairs <- list(
 # Liu et al. (2019): SI, CPD, AGESI, SC, DPW
  Liu_et_al_2019 = c("SI", "CPD", "AGESI", "SC", "DPW"),
 # Becker et al. (2021): CANNABIS
  Becker_et_al_2021 = c("CANNABIS")
) #list: substance_pairs

substance_dict <- subset_GWAS_pairs(GWAS_dict, substance_pairs)
cat("\n[substance_dict]\n")
 substance_dict

# -------------------------
# [Sociodemographics]
# -------------------------
sociodemo_pairs <- list(
 # Okbay et al. (2022): EA
 Okbay_et_al_2022  = c("EA"),  
 # Becker et al. (2021): RISK, FAMSAT, FINSAT, FRIENDSAT, WORKSAT
 Becker_et_al_2021 = c("RISK", "FAMSAT", "FINSAT", "FRIENDSAT", "WORKSAT"),
 # NEALE v2: TDI, RELIGATT
 NEALE_v2          = c("TDI", "RELIGATT")
) #list: sociodemo_pairs

sociodemo_dict <- subset_GWAS_pairs(GWAS_dict, sociodemo_pairs)
cat("\n[sociodemo_dict]\n")
 sociodemo_dict

# -------------------------
# [Psychiatric and Medical Conditions]
# -------------------------
medical_pairs <- list(
 # Becker et al. (2021): BMI
 Becker_et_al_2021 = c("BMI"),
 # NEALE v2: ASTHMA, COPD, HAYFEVER, HDL, BRON, DIABII, STROKE, HAMI, LUNG, HYPTENS, LONELY
 NEALE_v2          = c("ASTHMA", "COPD", "HAYFEVER", "HDL", "BRON", 
                       "DIABII", "STROKE", "HAMI", "LUNG", "HYPTENS", "LONELY"),
 # Gupta et al. (2024): NEURO, AGREE, CONSC, OPEN, EXTRA
 Gupta_et_al_2024  = c("NEURO", "AGREE", "CONSC", "OPEN", "EXTRA"),
 # PGC: SCHIZ, ADHD, MDD, ASD, BIP
 PGC               = c("SCHIZ", "ADHD", "MDD", "ASD", "BIP")
) #list: medical_pairs

medical_dict <- subset_GWAS_pairs(GWAS_dict, medical_pairs)
cat("\n[medical_dict]\n")
 medical_dict


# --------------------------------
# set trait list
# --------------------------------

for (TRAIT  in c("substance", "sociodemo", "medical")) {
    
    cat("\n*------------------- [",TRAIT,"] -------------------*\n")
    
    # get the corresponding dictionary
    dict <- get(paste0(TRAIT,"_dict"))
    
    # initialize trait path list
    trait_paths <- c()

    for (GWAS in names(dict)){

        phenos <- dict[[GWAS]]

        for (PHENO in phenos) {

            path <- paste0(org_dir, 
                          "/Other_Sumstats/",PHENO,"/",GWAS,
                          "/munge_sumstats/",PHENO,"_",GWAS,".sumstats.gz")

            trait_paths <- c(trait_paths, path)

        } #for: PHENO
    } #for: GWAS
    
    # assign back to each trait path list
    trait_path_name <- paste0("trait_",TRAIT)
    assign(trait_path_name, trait_paths)
     cat("\n[trait_",TRAIT,"]\n")
     print(get(trait_path_name))
    
    # assign trait names for each trait list
    trait_names <- paste0("trait_names_",TRAIT)
    assign(trait_names, unique(unlist(dict, use.names=F)))
     cat("\n[trait_names_",TRAIT,"]\n")
     print(get(trait_names))
    cat("\n*---------------------------------------------------*\n")
} #for: TRAIT

# birth cohort path
trait_birth  <-   paste0(root_dir, "/Meta_Analysis/",type,"/METAL_",ANCESTRY,"_",type,
                  "_MAresults_groupingBIN",1:4,".sumstats.gz")

# cobmine path
for (TRAIT  in c("substance", "sociodemo", "medical")) {
    
    trait_pheno <- get(paste0("trait_",TRAIT))
    assign(paste0("trait_input_",TRAIT,"_birth"), c(trait_pheno, trait_birth))
    
} #for: TRAIT

# --------------------------------
# Name the traits
# --------------------------------
trait_names_birth <- paste0("groupingBIN", 1:4)
for (TRAIT  in c("substance", "sociodemo", "medical")) {
    
    trait_pheno_names <- get(paste0("trait_names_",TRAIT))
    assign(paste0("trait_names_",TRAIT,"_birth"), c(trait_pheno_names, trait_names_birth))
    
} #for: TRAIT

# --------------------------------
# Define folder of LD scores
# --------------------------------
ld <- paste0(org_dir, "/", ANCESTRY, "/Source_Files/eur_w_ld_chr/")

# --------------------------------
# Define folder of LD weights [typically the same as folder of LD scores]
# --------------------------------
wld <- paste0(org_dir, "/", ANCESTRY, "/Source_Files/eur_w_ld_chr/")

# --------------------------------
# Define output directory
# --------------------------------
ldsc_out_dir <- paste0(root_dir, "/AR1_Model/results/", type,"/")
setwd(ldsc_out_dir)
getwd()

for (TRAIT  in c("substance", "sociodemo", "medical")) {
    cat("\n*------------------- [",TRAIT,"] -------------------*\n")

    # get input
    trait.input <- get(paste0("trait_input_",TRAIT,"_birth"))
    trait.name  <- get(paste0("trait_names_",TRAIT,"_birth"))
    
    # run LDSC
    LDSCoutput <- ldsc(traits = trait.input,
                       sample.prev = NA,
                       population.prev = NA,
                       ld = ld,
                       wld = wld,
                       trait.names = trait.name
                      ) #ldsc
    
    # assign results
    output_name <- paste0("LDSCoutput_",TRAIT)
    assign(output_name, LDSCoutput)
    
    cat("\n*---------------------------------------------------*\n")   
} #for: TRAIT

for (TRAIT  in c("substance", "sociodemo", "medical")) {
    cat("\n*------------------- [",TRAIT,"] -------------------*\n")

    # get LDSCoutput
    LDSCoutput <- get(paste0("LDSCoutput_",TRAIT))
    
    # Plot genetic covariance matrix (S)
    corrplot::corrplot(cov2cor(LDSCoutput$S), 
                       type = "upper",
                       tl.col = "black",
                       tl.cex = 0.8,
                       is.corr = F)
    #method = "color",
    cat("\n*---------------------------------------------------*\n")   
} #for: TRAIT

# A =~ B: Build latent variable; A (latent factor) is measured by B (observed variable)
# A ~  B: Regression; Regression of A on B
# A ~~ B: Correlation; Genetic correlation or covariance between A and B

build_LCSM <- function(n_start = 1, n_end = 4, external_traits = NULL) {
  
  # sequences
  bins <- n_start:n_end
  bins_prev   <- n_start:(n_end-1)
  bins_curr   <- (n_start+1):n_end
  
  # -------------------------
  # Measurement layer
  # - binX is an observed GWAS trait
  # - BINX is a latent factor representing its genetic component; total genetic signal at time t
  # -------------------------
  meas <- paste0(
    "BIN", bins, " =~ 1*groupingBIN", bins,
    collapse = "\n"
  )
  
  # -------------------------
  # Innovation (Genetic change) layer
  # - DX is a latent innovation factor; new genetic signal at time t
  # - D5,..., D19 are innovations meaning the new genetic variance at that wave after removing what is carried over from the previous wave.
  # -------------------------
  D_layer <- paste0(
    "D", bins_curr, " =~ 1*BIN", bins_curr,
    collapse = "\n"
  )
  
  # -------------------------
  # Autoregressive structure
  # -------------------------
  AR <- paste0(
    "BIN", bins_curr, " ~ BIN", bins_prev,
    collapse = "\n"
  )
  
  # -------------------------
  # D variances
  # -------------------------
  D_var <- paste0(
    "D", bins_curr, " ~~ vD", bins_curr, "*D", bins_curr,
    collapse = "\n"
  )
  
      # -------------------------
      # Positivity constraints on all D variances
      # -------------------------
      D_pos <- paste0(
        "vD", bins_curr, " > 0.0000001",
        collapse = "\n"
      )
  
  # -------------------------
  #  Innovation independence
  # -------------------------
  D_indep <- ""
  for (i in bins_curr[-length(bins_curr)]) {
    rhs <- paste0("0*D", (i+1):n_end, collapse = " + ")
    D_indep <- paste0(D_indep, "D", i, " ~~ ", rhs, "\n")
  }
  
  # -------------------------
  # latent variables independent conditional on AR structure
  # -------------------------
  BIN_indep <- ""
  for (i in bins[-length(bins)]) {
    rhs <- paste0("0*BIN", (i+1):n_end, collapse = " + ")
    BIN_indep <- paste0(BIN_indep, "BIN", i, " ~~ ", rhs, "\n")
  }
  
  # -------------------------
  # Baseline orthogonality
  # -------------------------
  baseline_ortho <- paste0(
    "BIN", n_start, " ~~ ",
    paste0("0*D", bins_curr, collapse = " + ")
  )
  
  # -------------------------
  # External traits
  # -------------------------
  #Interpretation example: ea ~~ D10 = 0.25 means genetic variants associated with EA are positively associated with increase in smoking initiation in bin10
  #                                  =-0.30 means EA genetics are protective. That is, they reduce increase in smoking at that time.
  external_block <- ""
  
  if (!is.null(external_traits)) {
    
    D_terms <- paste0("D", bins_curr, collapse = " + ")
    
    for (trait in external_traits) {
      
      # baseline
      external_block <- paste0(
        external_block,
        trait, " ~~ BIN", n_start, "\n"
      )
      
      # time-specific innovation associations
      external_block <- paste0(
        external_block,
        trait, " ~~ ", D_terms, "\n"
      )
    }
  }
  
  # -------------------------
  # Combine everything
  # -------------------------
  model <- paste0(
"
# -------------------------
# Measurement layer
# -------------------------
", meas, "

# -------------------------
# Innovation layer
# -------------------------
", D_layer, "

# -------------------------
# Autoregressive structure
# -------------------------
", AR, "

# -------------------------
# D variances
# -------------------------
", D_var, "

# Positivity constraints
", D_pos, "

# -------------------------
# Innovations independent
# -------------------------
", D_indep, "

# -------------------------
# Latent variable independence
# -------------------------
", BIN_indep, "

# -------------------------
# Baseline orthogonality
# -------------------------
", baseline_ortho, "

# -------------------------
# External traits
# -------------------------
", external_block, "
"
  )
  
  return(model)
}

for (TRAIT  in c("substance", "sociodemo", "medical")) {
    cat("\n*------------------- [",TRAIT,"] -------------------*\n")

    # define names
    trait_names <- paste0("trait_names_", TRAIT)
    
    # build model
    LCSM <- build_LCSM(
                                 n_start = 1,
                                 n_end = 4,
                                 external_traits = get(trait_names)
                                 ) # build LCSM
    
    # assign model
    model_name  <- paste0("LCSM_", TRAIT)
    assign(model_name, LCSM)
    cat(get(model_name))

    cat("\n*---------------------------------------------------*\n")   
} #for: TRAIT

for (TRAIT  in c("substance", "sociodemo", "medical")) {
    cat("\n*------------------- [",TRAIT,"] -------------------*\n")

    # define names
    LDSCoutput_name <- paste0("LDSCoutput_", TRAIT)
    LCSM_model_name <- paste0("LCSM_", TRAIT)
    
    # fit the structural equation model using the LDSC covariance matrix with GenomicSEM.
    LCSM.fit<- usermodel(covstruc = get(LDSCoutput_name),
                            estimation = "DWLS",
                            model = get(LCSM_model_name),
                            imp_cov = T)

    # assign results
    LCSM_fit_name <- paste0("LCSM.fit.",TRAIT)
    assign(LCSM_fit_name, LCSM.fit)
    
    cat("\n[LCSM.fit.",TRAIT,"$modelfit]\n")
    print(get(LCSM_fit_name)$modelfit)
    cat("\n[LCSM.fit.",TRAIT,"$results]\n")
    print(get(LCSM_fit_name)$results)
    # Unstand_Est (unstandardized estimate): Raw covariance estimate; Genetic covariance
    # Unstand_SE (standard error): Uncertainty of the estimate
    # STD_Genotype: Genetic correlation; Standardized genetic covariance
    # STD_Genotype_SE: Standard error of the genetic correlation
    # STD_All: Fully standardized estimate; Standardizes both latent variables and observed variables
    # p_value: significance
    
    cat("\n*---------------------------------------------------*\n")   
} #for: TRAIT

for (TRAIT  in c("substance", "sociodemo", "medical")) {
    cat("\n*------------------- [",TRAIT,"] -------------------*\n")

    # define names
    results_name <- paste0("LCSM.fit.", TRAIT)
    trait_name   <- paste0("trait_names_", TRAIT)
    
    results <- get(results_name)$results
     #head(results, n=50)
    
    # filter results
    results_filtered <- results %>%
                            filter(
                                op == "~~",                   # correlations
                                rhs %in% get(trait_name),     # trait on LHS
                                grepl("^D[0-9]+$", lhs)       # D5-D19 on RHS
                            )
    
    # assign
    results_filtered_name <- paste0("results_",TRAIT)
    assign(results_filtered_name, results_filtered)
     cat("\n[results_",TRAIT,"]\n")
     print(get(results_filtered_name))
    
    cat("\n*---------------------------------------------------*\n")   
} #for: TRAIT

# --------------------------------
# Dictionary (abb -> full name)
# --------------------------------
#define [pheno_name_dict] with [Names=Abbreviation] and [Values=full name]
pheno_name_dict <- list(
                        EA = "Educational Attainment",
                        NEURO = "Big 5 Personality: Neuroticism",
                        AGREE = "Big 5 Personality: Agreeableness", 
                        CONSC = "Big 5 Personality: Conscientiousness",
                        OPEN = "Big 5 Personality: Openness", 
                        EXTRA = "Big 5 Personality: Extraversion",
                        RISK = "Risk Tolerance",
                        FAMSAT = "Life Satisfaction: Family", 
                        FINSAT = "Life Satisfaction: Finance", 
                        FRIENDSAT = "Life Satisfaction: Friends", 
                        WORKSAT = "Life Satisfaction: Work",
                        LONELY = "Loneliness",
                        RELIGATT = "Religious Attendance",
                        SI = "Smoking Initiation",
                        CPD = "Cigarettes per Day",
                        AGESI = "Age of Smoking Initiation",
                        SC = "Smoking Cessation",
                        DPW = "Drinks per Week",
                        CANNABIS = "Cannabis Use",
                        TDI = "Townsend Deprivation Index",
                        BMI = "Body Mass Index",
                        ASTHMA = "Asthma",
                        COPD = "Chronic Obstructive Pulmonary Disease",
                        HAYFEVER = "Hayfever/Allergic Rhinitis",
                        HDL = "High Cholesterol",
                        DIABII = "Type II Diabetes",
                        STROKE = "Stroke",
                        HAMI = "Heart Attack/Myocardial Infarction",
                        LUNG = "Lung Cancer",
                        HYPTENS = "Hypertension",
                        SCHIZ = "Schizophrenia",
                        ADHD = "Attention Deficit Hyperactivity Disorder",
                        MDD = "Major Depressive Disorder",
                        ASD = "Autism Spectrum Disorder",
                        BIP = "Bipolar Disorder"
                       ) #list: pheno_name_dict (Names=GWAS, Values=phenotypes)

# convert list to named character vector
pheno_vec <- unlist(pheno_name_dict)

# Generate [Phenotype column]
## Substance Use
results_substance$Phenotype <- ifelse(
                                      results_substance$rhs %in% names(pheno_vec),
                                      pheno_vec[results_substance$rhs],
                                      results_substance$rhs
                                     ) #ifelse

## Sociodemographic
results_sociodemo$Phenotype <- ifelse(
                                      results_sociodemo$rhs %in% names(pheno_vec),
                                      pheno_vec[results_sociodemo$rhs],
                                      results_sociodemo$rhs
                                     ) #ifelse

## Psychiatric and medical conditions
results_medical$Phenotype <- ifelse(
                                      results_medical$rhs %in% names(pheno_vec),
                                      pheno_vec[results_medical$rhs],
                                      results_medical$rhs
                                     ) #ifelse

# combine data
data <- rbind(results_substance, results_sociodemo, results_medical)
data

# check min and max of [STD_Genotype_SE]
cat("\nMin of [STD_Genotype_SE]\n")
print(min(data$STD_Genotype))
cat("\nMax of [STD_Genotype_SE]\n")
print(max(data$STD_Genotype))

n_D <- length(unique(data$lhs))
n_D

# -------------------------
# Block labels
# -------------------------
data$Block <- NA
data$Block[data$rhs %in% trait_names_substance] <- "Substance Use"
data$Block[data$rhs %in% trait_names_sociodemo] <- "Sociodemographics"
data$Block[data$rhs %in% trait_names_medical]   <- "Psychiatric and Medical Conditions"

# make facet order explicit
data$Block <- factor(data$Block, levels = c("Substance Use","Sociodemographics","Psychiatric and Medical Conditions"))
# make innovation order explicit
data$lhs <- trimws(as.character(data$lhs))
# make innovation order explicit
data$lhs <- trimws(as.character(data$lhs))

data$lhs <- factor(
  data$lhs,
  levels = c("D2", "D3", "D4"),
  labels = c("Baby Boomer-specific",
             "Generation X-specific",
             "Millennial-specific")
)# factor

# make phenotye order explicit
## Preserve original order
pheno_order <- c("Smoking Initiation", "Cigarettes per Day", "Age of Smoking Initiation",
                 "Smoking Cessation", "Drinks per Week", "Cannabis Use",
                 "Educational Attainment", "Risk Tolerance", "Life Satisfaction: Family",
                 "Life Satisfaction: Finance", "Life Satisfaction: Friends", "Life Satisfaction: Work",
                 "Townsend Deprivation Index", "Religious Attendance",
                 "Body Mass Index", "Asthma", "Chronic Obstructive Pulmonary Disease", 
                 "Hayfever/Allergic Rhinitis", "High Cholesterol",
                 "Type II Diabetes", "Stroke", "Heart Attack/Myocardial Infarction", 
                 "Lung Cancer", "Hypertension", "Loneliness",
                 "Big 5 Personality: Neuroticism", "Big 5 Personality: Agreeableness",
                 "Big 5 Personality: Conscientiousness", "Big 5 Personality: Openness",
                 "Big 5 Personality: Extraversion",
                 "Schizophrenia", "Attention Deficit Hyperactivity Disorder", "Major Depressive Disorder",
                 "Autism Spectrum Disorder", "Bipolar Disorder"
                 )
pheno_order <- rev(pheno_order) # top-to-bottom order
data$Phenotype <- factor(data$Phenotype, levels = pheno_order)

# -------------------------
# Test for Significance ([rg=0]???)
#  with Group specific bonferroni correction.
# -------------------------
cat("\nNumber of Tests: \n")
n_pheno <- length(unique(data$Phenotype))
n_pheno
cat("\nNumber of Innovations: \n")
n_D <- length(unique(data$lhs))
n_D

data <- data %>%
          mutate(
            num_tests = n_D*n_pheno,
            #z = rg / SE,
            #p = 2 * pnorm(-abs(z)),
            p_value = as.numeric(p_value),
            stars = case_when(
              p_value < (0.05 / num_tests) ~ "**",
              p_value < 0.05 ~ "*",
              TRUE ~ ""
            ) #case_when
          ) #mutate

#data
# ------------------------
# 5) Plot
# -------------------------
pheno_inno <- ggplot(data, aes(lhs, Phenotype, fill = STD_Genotype)) +
              geom_tile(color = "white", linewidth = 0.7) +
              geom_text(aes(label = paste0(sprintf("%.2f", STD_Genotype), stars)), size = 2.9) +
              scale_fill_gradient2(
              low = "#A020F0",   # bright purple
              mid = "white",
              high = "#FF0000",  # bright red
                midpoint = 0,
                limits = c(-0.25, 0.25),
                name = "Genetic Correlation",
                na.value = "white"
              ) +
              facet_grid(Block ~ ., scales = "free_y", space = "free_y") +
              scale_x_discrete(position = "top") +
              scale_y_discrete(position = "right") +
              theme_minimal(base_size = 12) +
              theme(
                panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text.x = element_text(size = 8),
                axis.text.y.right = element_text(hjust = 0),
                strip.background.y = element_rect(fill = "grey90", color = "black"),
                strip.text.y = element_text(angle = 90, face = "bold", size = 7),
                legend.position = "none"
              )
pheno_inno

# Save it
ggsave(paste0(output_dir, "/",type,"/Figures/pheno_inno_all_gen.png"), plot=pheno_inno, width = 10, height = 6, dpi = 500)


