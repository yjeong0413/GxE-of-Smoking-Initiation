R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()

cat("\nLoad Libraries\n")
cat("\n-----------------------------------------------------------------------------------\n")
cat("\n*Load [data.table]\n")
library(data.table)
cat("\n*Load [dplyr]\n")
library(dplyr)
cat("\n")
cat("\n*Load [stringr]\n")
library(stringr)
cat("\n")
cat("\n*Load [tidyr]\n")
library(tidyr)
cat("\n")
#Libraries useful for creating latex tables (.tex file)
cat("\n*Load [stargazer]\n") 
 #install.packages("stargazer")
library(stargazer)
cat("\n")
cat("\n*Load [kableExtra]\n") 
 #install.packages("kableExtra")
#library(kableExtra)
cat("\n")
cat("\n*Load [xtable]\n") 
library(xtable)
cat("\n")
#cat("\n*Load [knitr]\n") 
#library(knitr)
cat("\n")
#cat("\n*Load [tidyverse]\n") 
#library(tidyverse)
cat("\n")
cat("\n*Load [ggplot2]\n") 
library(ggplot2)
cat("\n")
cat("\n*Load [reshape2]\n") 
#install.packages("reshape2")
library(reshape2)
cat("\n")
cat("\n*Load [purrr]\n") 
library(purrr)
cat("\n")
cat("\n*Load [tibble]\n") 
library(tibble)
cat("\n")
cat("\n*Load [patchwork]\n") 
library(patchwork)
cat("\n")
cat("\n*Load [gt]\n") 
library(gt)
cat("\n")
cat("\n*Load [writexl]\n") 
library(writexl)
cat("\n")
cat("\n-----------------------------------------------------------------------------------\n")


#Basic settings
cat("\nBasic setting\n")
ANCESTRY="european"
cat("\n-------------------------------------------------------------------------------\n")

#Define directories
cat("\nSet directories\n")
root_dir=paste0("/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc/", ANCESTRY ,"/CLEAN/GWAS_QC")
working_dir=paste0(root_dir,"/Genetic_Correlation/raw")

cat("\n-------------------------------------------------------------------------------\n")


# --------------------------------
# Sample
# --------------------------------
sample="Full"

# --------------------------------
# Dictionary (GWAS -> phenotypes)
# --------------------------------
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
# Storage object (clean structure)
# --------------------------------
all_results <- list()
all_results_se <- list()

# --------------------------------
# Loop over GWAS
# --------------------------------
for (GWAS in names(GWAS_dict)) {
    
    PHENOTYPE_LIST <- GWAS_dict[[GWAS]]
    all_matrices <- list()
    all_matrices_se <- list()
    
    cat("\n# ------------------ ", GWAS, " ------------------# \n")
    
    # --------------------------------
    # Loop over phenotypes
    # --------------------------------    
    for (PHENOTYPE in PHENOTYPE_LIST) {
    # Set working directory
    results_dir=paste0(working_dir,"/", PHENOTYPE,"/", GWAS)
    results_dir_tables=paste0(results_dir,"/Tables")
    results_dir_figures=paste0(results_dir,"/Figures")

    #print(results_dir)
    # ---------------- Setup ----------------
    dir.exists(file.path(results_dir))
    setwd(file.path(results_dir, "/full"))  # Change to where your .log files are stored

    # List all LDSC log files matching pattern
    log_files <- list.files(pattern = "_rg_raw.log$")

    # ---------------- Extraction Functions ----------------
    # Extraction function for rg summary
    extract_rg_summary <- function(file) {
      lines <- readLines(file)
      start <- grep("Summary of Genetic Correlation Results", lines)
      if (length(start) == 0) return(NULL)
  
      data_lines <- lines[(start + 1):length(lines)]
      data_lines <- data_lines[data_lines != ""]

      tryCatch({
        df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
        colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                          "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                          "gcov_int", "gcov_int_se")
    
        # Labeling function for phenotype identity (customize this per use case)
        extract_label <- function(path) {
          if (grepl(GWAS, path)) return(PHENOTYPE)
          else if (grepl("ever_tobacco_user", path)) return("Smoking")
          else return("Other")
        }
    
        df$Trait1 <- sapply(df$p1, extract_label)
        df$Trait2 <- sapply(df$p2, extract_label)

        return(df)
      }, error = function(e) NULL)
    }

    # ---------------- Combine All Summaries ----------------
    all_data <- map_df(log_files, extract_rg_summary)

    # ---------------- Construct Matrix ----------------
    traits <- c(PHENOTYPE, "Smoking")
    ## rg
    rg_matrix <- all_data %>%
                  filter(Trait1 %in% traits, Trait2 %in% traits) %>%
                  distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                  select(Trait1, Trait2, rg) %>%
                  tidyr::pivot_wider(names_from = Trait2, values_from = rg) %>%
                  tibble::column_to_rownames("Trait1") %>%
                  dplyr::select(any_of(traits)) %>%
                  .[intersect(rownames(.), traits), , drop = FALSE]
    ## rg_se
    rg_se_matrix <- all_data %>%
                  filter(Trait1 %in% traits, Trait2 %in% traits) %>%
                  distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                  select(Trait1, Trait2, se) %>%
                  tidyr::pivot_wider(names_from = Trait2, values_from = se) %>%
                  tibble::column_to_rownames("Trait1") %>%
                  dplyr::select(any_of(traits)) %>%
                  .[intersect(rownames(.), traits), , drop = FALSE]

    # ---------------- Output ----------------
    rownames(rg_matrix)    <- make.names(rownames(rg_matrix))
    rownames(rg_se_matrix) <- make.names(rownames(rg_se_matrix))
    colnames(rg_matrix)    <- make.names(colnames(rg_matrix))
    colnames(rg_se_matrix) <- make.names(colnames(rg_se_matrix))

    # Print matrix
    #print(rg_matrix)

    # Save results in a list: all_matrices
    rg_matrix_df <- as.data.frame(rg_matrix)
    rg_se_matrix_df <- as.data.frame(rg_se_matrix)
    #rg_matrix_df$Phenotype <- PHENOTYPE # keep phenotype label
    all_matrices[[PHENOTYPE]] <- rg_matrix_df
    all_matrices_se[[PHENOTYPE]] <- rg_se_matrix_df

    } #for: PHENOTYPE

    # --------------------------------
    # Combine phenotypes for this GWAS
    # --------------------------------
    # rg
    cat("\nFinal Genetic Correaltion Matrix:\n")
    mtx_name <- paste("GC_", sample,"_", GWAS)
    final_df <- dplyr::bind_rows(all_matrices)
    colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name, final_df)
    print(get(mtx_name))
    # rg_se
    cat("\nFinal SE of Genetic Correaltion Matrix:\n")
    mtx_name_se <- paste("GC_se_", sample,"_", GWAS)
    final_df_se <- dplyr::bind_rows(all_matrices_se)
    colnames(final_df_se)[colnames(final_df_se) == "Smoking"] <- sample
    assign(mtx_name_se, final_df_se)
    print(get(mtx_name_se))
    
    # Store directly
    ## rg
    all_results[[sample]][[GWAS]] <- get(mtx_name)
    ## rg_se
    all_results_se[[sample]][[GWAS]] <- get(mtx_name_se)
           
} #for: GWAS

# --------------------------------
# Sample
# --------------------------------
sample="gender"

# --------------------------------
# Dictionary (GWAS -> phenotypes)
# --------------------------------
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
# Storage object (clean structure)
# --------------------------------
#all_results <- list()

# --------------------------------
# Loop over GWAS
# --------------------------------
for (GWAS in names(GWAS_dict)) {
    
    PHENOTYPE_LIST <- GWAS_dict[[GWAS]]
    all_matrices <- list()
    all_matrices_se <- list()
    
    cat("\n# ------------------ ", GWAS, " ------------------# \n")
    
    # --------------------------------
    # Loop over phenotypes
    # --------------------------------    
    for (PHENOTYPE in PHENOTYPE_LIST) {
    # Set working directory
    results_dir=paste0(working_dir,"/", PHENOTYPE,"/", GWAS)
    results_dir_tables=paste0(results_dir,"/Tables")
    results_dir_figures=paste0(results_dir,"/Figures")

    #print(results_dir)
    # ---------------- Setup ----------------
    dir.exists(file.path(results_dir))
    setwd(file.path(results_dir, "/gender"))  # Change to where your .log files are stored

    # List all LDSC log files matching pattern
    log_files <- list.files(pattern = "_rg_raw.log$")

    # ---------------- Extraction Functions ----------------
    # Extraction function for rg summary
    extract_rg_summary <- function(file) {
      lines <- readLines(file)
      start <- grep("Summary of Genetic Correlation Results", lines)
      if (length(start) == 0) return(NULL)
  
      data_lines <- lines[(start + 1):length(lines)]
      data_lines <- data_lines[data_lines != ""]

      tryCatch({
        df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
        colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                          "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                          "gcov_int", "gcov_int_se")
          
        # Custom label extraction for EA vs gender
        extract_gender <- function(path) {
          if (grepl("female", path)) return("female")
          else if (grepl("male", path)) return("male")
          else if (grepl(PHENOTYPE, path)) return(PHENOTYPE)
          else return(NA)
        }
        df$Trait1 <- sapply(df$p1, extract_gender)
        df$Trait2 <- sapply(df$p2, extract_gender)

        return(df)
      }, error = function(e) NULL)
    }

    # ---------------- Combine All Summaries ----------------
    all_data <- map_df(log_files, extract_rg_summary)

    # ---------------- Construct Matrix ----------------
    traits <- c(PHENOTYPE, "female", "male")
    # rg
    rg_matrix_gender <- all_data %>%
                          filter(Trait1 %in% traits, Trait2 %in% traits) %>%
                          distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                          select(Trait1, Trait2, rg) %>%
                          pivot_wider(names_from = Trait2, values_from = rg) %>%
                          column_to_rownames("Trait1") %>%
                          select(any_of(traits)) %>%
                          .[intersect(rownames(.), traits), , drop = FALSE]
    # rg
    rg_se_matrix_gender <- all_data %>%
                          filter(Trait1 %in% traits, Trait2 %in% traits) %>%
                          distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                          select(Trait1, Trait2, se) %>%
                          pivot_wider(names_from = Trait2, values_from = se) %>%
                          column_to_rownames("Trait1") %>%
                          select(any_of(traits)) %>%
                          .[intersect(rownames(.), traits), , drop = FALSE]

    # ---------------- Output ----------------
    rownames(rg_matrix_gender)    <- make.names(rownames(rg_matrix_gender))
    rownames(rg_se_matrix_gender) <- make.names(rownames(rg_se_matrix_gender))
    colnames(rg_matrix_gender)    <- make.names(colnames(rg_matrix_gender))
    colnames(rg_se_matrix_gender) <- make.names(colnames(rg_se_matrix_gender))

    # Print matrix
    #print(rg_matrix)

    # Save results in a list: all_matrices
    rg_matrix_df <- as.data.frame(rg_matrix_gender)
    rg_se_matrix_df <- as.data.frame(rg_se_matrix_gender)
    #rg_matrix_df$Phenotype <- PHENOTYPE # keep phenotype label
    all_matrices[[PHENOTYPE]] <- rg_matrix_df
    all_matrices_se[[PHENOTYPE]] <- rg_se_matrix_df

    } #for: PHENOTYPE

    # --------------------------------
    # Combine phenotypes for this GWAS
    # --------------------------------
    # rg
    cat("\nFinal Genetic Correaltion Matrix:\n")
    mtx_name <- paste("GC_", sample,"_", GWAS)
    final_df <- dplyr::bind_rows(all_matrices)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name, final_df)
    print(get(mtx_name))
    # rg_se
    cat("\nFinal SE of Genetic Correaltion Matrix:\n")
    mtx_name_se <- paste("GC_se_", sample,"_", GWAS)
    final_df_se <- dplyr::bind_rows(all_matrices_se)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name_se, final_df_se)
    print(get(mtx_name_se))
    
    # Store directly
    # rg
    all_results[[sample]][[GWAS]] <- get(mtx_name)
    # rg_se
    all_results_se[[sample]][[GWAS]] <- get(mtx_name_se)
           
} #for: GWAS

# --------------------------------
# Sample
# --------------------------------
sample="region"

# --------------------------------
# Dictionary (GWAS -> phenotypes)
# --------------------------------
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
# Storage object (clean structure)
# --------------------------------
#all_results <- list()

# --------------------------------
# Loop over GWAS
# --------------------------------
for (GWAS in names(GWAS_dict)) {
    
    PHENOTYPE_LIST <- GWAS_dict[[GWAS]]
    all_matrices <- list()
    all_matrices_se <- list()
    
    cat("\n# ------------------ ", GWAS, " ------------------# \n")
    
    # --------------------------------
    # Loop over phenotypes
    # --------------------------------    
    for (PHENOTYPE in PHENOTYPE_LIST) {
    # Set working directory
    results_dir=paste0(working_dir,"/", PHENOTYPE,"/", GWAS)
    results_dir_tables=paste0(results_dir,"/Tables")
    results_dir_figures=paste0(results_dir,"/Figures")

    #print(results_dir)
    # ---------------- Setup ----------------
    dir.exists(file.path(results_dir))
    setwd(file.path(results_dir, "/region"))  # Change to where your .log files are stored

    # List all LDSC log files matching pattern
    log_files <- list.files(pattern = "_rg_raw.log$")

    # ---------------- Extraction Functions ----------------
    # Function to extract region label from file path
    extract_region <- function(path) {
      if (grepl("MW/", path)) return("MW")
      else if (grepl("NE/", path)) return("NE")
      else if (grepl("SE/", path)) return("SE")
      else if (grepl("SW/", path)) return("SW")
      else if (grepl("W/", path)) return("W")
      else if (grepl(PHENOTYPE, path) || grepl(GWAS, path)) return(PHENOTYPE)
      else return(NA)
    }
    # Extraction function for rg summary
    extract_rg_summary <- function(file) {
      lines <- readLines(file)
      start <- grep("Summary of Genetic Correlation Results", lines)
      if (length(start) == 0) return(NULL)
  
      data_lines <- lines[(start + 1):length(lines)]
      data_lines <- data_lines[data_lines != ""]

      tryCatch({
        df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
        colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                          "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                          "gcov_int", "gcov_int_se")          
        
        df$Trait1 <- sapply(df$p1, extract_region)
        df$Trait2 <- sapply(df$p2, extract_region)

        return(df)
      }, error = function(e) NULL)
    }

    # ---------------- Combine All Summaries ----------------
    all_data <- map_df(log_files, extract_rg_summary)

    # ---------------- Construct Matrix ----------------
    # (Optional) If you're comparing across multiple traits
     regions <- c("MW", "NE", "SE", "SW", "W")
    # rg
    rg_matrix_region <- all_data %>%
                          filter(Trait1 == PHENOTYPE & Trait2 %in% regions) %>%
                          distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                          select(Trait1, Trait2, rg) %>%
                          pivot_wider(names_from = Trait2, values_from = rg) %>%
                          column_to_rownames("Trait1")
    # rg_se
    rg_se_matrix_region <- all_data %>%
                          filter(Trait1 == PHENOTYPE & Trait2 %in% regions) %>%
                          distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                          select(Trait1, Trait2, se) %>%
                          pivot_wider(names_from = Trait2, values_from = se) %>%
                          column_to_rownames("Trait1")

    # ---------------- Output ----------------
    # Optional: make column names safe
    colnames(rg_matrix_region) <- make.names(colnames(rg_matrix_region))
    colnames(rg_se_matrix_region) <- make.names(colnames(rg_se_matrix_region))

    # Print matrix
    #print(rg_matrix)

    # Save results in a list: all_matrices
    rg_matrix_df <- as.data.frame(rg_matrix_region)
    rg_se_matrix_df <- as.data.frame(rg_se_matrix_region)
    #rg_matrix_df$Phenotype <- PHENOTYPE # keep phenotype label
    all_matrices[[PHENOTYPE]] <- rg_matrix_df
    all_matrices_se[[PHENOTYPE]] <- rg_se_matrix_df

    } #for: PHENOTYPE

    # --------------------------------
    # Combine phenotypes for this GWAS
    # --------------------------------
    # rg
    cat("\nFinal Genetic Correaltion Matrix:\n")
    mtx_name <- paste("GC_", sample,"_", GWAS)
    final_df <- dplyr::bind_rows(all_matrices)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name, final_df)
    print(get(mtx_name))
    # rg_se
    cat("\nFinal SE of Genetic Correaltion Matrix:\n")
    mtx_name_se <- paste("GC_se_", sample,"_", GWAS)
    final_df_se <- dplyr::bind_rows(all_matrices_se)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name_se, final_df_se)
    print(get(mtx_name_se))
    
    # Store directly
    # rg
    all_results[[sample]][[GWAS]] <- get(mtx_name)
    # rg_se
    all_results_se[[sample]][[GWAS]] <- get(mtx_name_se)
           
} #for: GWAS

# --------------------------------
# Sample
# --------------------------------
sample="birthyear"

# --------------------------------
# Dictionary (GWAS -> phenotypes)
# --------------------------------
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
# Storage object (clean structure)
# --------------------------------
#all_results <- list()

# --------------------------------
# Loop over GWAS
# --------------------------------
for (GWAS in names(GWAS_dict)) {
    
    PHENOTYPE_LIST <- GWAS_dict[[GWAS]]
    all_matrices <- list()
    all_matrices_se <- list()
    
    cat("\n# ------------------ ", GWAS, " ------------------# \n")
    
    # --------------------------------
    # Loop over phenotypes
    # --------------------------------    
    for (PHENOTYPE in PHENOTYPE_LIST) {
    # Set working directory
    results_dir=paste0(working_dir,"/", PHENOTYPE,"/", GWAS)
    results_dir_tables=paste0(results_dir,"/Tables")
    results_dir_figures=paste0(results_dir,"/Figures")

    #print(results_dir)
    # ---------------- Setup ----------------
    dir.exists(file.path(results_dir))
    setwd(file.path(results_dir, "/birthyear"))  # Change to where your .log files are stored

    # List all LDSC log files matching pattern
    log_files <- list.files(pattern = "_rg_raw.log$")

    # ---------------- Extraction Functions ----------------
    # Extraction function for rg summary
    extract_rg_summary <- function(file) {
      lines <- readLines(file)
      start <- grep("Summary of Genetic Correlation Results", lines)
      if (length(start) == 0) return(NULL)
  
      data_lines <- lines[(start + 1):length(lines)]
      data_lines <- data_lines[data_lines != ""]

      tryCatch({
        df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
        colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                          "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                          "gcov_int", "gcov_int_se")          
        
        df$Trait1 <- PHENOTYPE
        df$Trait2 <- sub(".*birthyear_bin(\\d+).*", "bin\\1", df$p2)

        return(df)
      }, error = function(e) NULL)
    }

    # ---------------- Combine All Summaries ----------------
    all_data <- map_df(log_files, extract_rg_summary)

    # ---------------- Construct Matrix ----------------
    # (Optional) If you're comparing across multiple traits
    bins <- paste0("bin",4:19)
    # rg
    rg_matrix_birthyear <- all_data %>%
                              filter(Trait1 == PHENOTYPE & Trait2 %in% bins) %>%
                              distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                              select(Trait1, Trait2, rg) %>%
                              pivot_wider(names_from = Trait2, values_from = rg) %>%
                              column_to_rownames("Trait1")
    # rg_se
    rg_se_matrix_birthyear <- all_data %>%
                              filter(Trait1 == PHENOTYPE & Trait2 %in% bins) %>%
                              distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                              select(Trait1, Trait2, se) %>%
                              pivot_wider(names_from = Trait2, values_from = se) %>%
                              column_to_rownames("Trait1")

    # ---------------- Output ----------------
    # Print matrix
    #print(rg_matrix)

    # Save results in a list: all_matrices
    rg_matrix_df <- as.data.frame(rg_matrix_birthyear)
    rg_se_matrix_df <- as.data.frame(rg_se_matrix_birthyear)
    #rg_matrix_df$Phenotype <- PHENOTYPE # keep phenotype label
    all_matrices[[PHENOTYPE]] <- rg_matrix_df
    all_matrices_se[[PHENOTYPE]] <- rg_se_matrix_df

    } #for: PHENOTYPE

    # --------------------------------
    # Combine phenotypes for this GWAS
    # --------------------------------
    # rg
    cat("\nFinal Genetic Correaltion Matrix:\n")
    mtx_name <- paste("GC_", sample,"_", GWAS)
    final_df <- dplyr::bind_rows(all_matrices)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name, final_df)
    print(get(mtx_name))
    # rg_se
    cat("\nFinal SE of Genetic Correaltion Matrix:\n")
    mtx_name_se <- paste("GC_se_", sample,"_", GWAS)
    final_df_se <- dplyr::bind_rows(all_matrices_se)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name_se, final_df_se)
    print(get(mtx_name_se))
    
    # Store directly
    # rg
    all_results[[sample]][[GWAS]] <- get(mtx_name)
    # rg_se
    all_results_se[[sample]][[GWAS]] <- get(mtx_name_se)
           
} #for: GWAS

# --------------------------------
# Sample
# --------------------------------
sample="generation"

# --------------------------------
# Dictionary (GWAS -> phenotypes)
# --------------------------------
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
# Storage object (clean structure)
# --------------------------------
#all_results <- list()

# --------------------------------
# Loop over GWAS
# --------------------------------
for (GWAS in names(GWAS_dict)) {
    
    PHENOTYPE_LIST <- GWAS_dict[[GWAS]]
    all_matrices <- list()
    all_matrices_se <- list()
    
    cat("\n# ------------------ ", GWAS, " ------------------# \n")
    
    # --------------------------------
    # Loop over phenotypes
    # --------------------------------    
    for (PHENOTYPE in PHENOTYPE_LIST) {
    # Set working directory
    results_dir=paste0(working_dir,"/", PHENOTYPE,"/", GWAS)
    results_dir_tables=paste0(results_dir,"/Tables")
    results_dir_figures=paste0(results_dir,"/Figures")

    #print(results_dir)
    # ---------------- Setup ----------------
    dir.exists(file.path(results_dir))
    setwd(file.path(results_dir, "/birthyear"))  # Change to where your .log files are stored

    # List all LDSC log files matching pattern
    log_files <- list.files(pattern = "_rg_raw_2ndMA.log$")

    # ---------------- Extraction Functions ----------------
    # Extraction function for rg summary
    extract_rg_summary <- function(file) {
      lines <- readLines(file)
      start <- grep("Summary of Genetic Correlation Results", lines)
      if (length(start) == 0) return(NULL)
  
      data_lines <- lines[(start + 1):length(lines)]
      data_lines <- data_lines[data_lines != ""]

      tryCatch({
        df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
        colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                          "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                          "gcov_int", "gcov_int_se")          
        
        df$Trait1 <- PHENOTYPE
        df$Trait2 <- sub(".*_groupingBIN(\\d+).*", "groupingBIN\\1", df$p2)

        return(df)
      }, error = function(e) NULL)
    }

    # ---------------- Combine All Summaries ----------------
    all_data <- map_df(log_files, extract_rg_summary)

    # ---------------- Construct Matrix ----------------
    # (Optional) If you're comparing across multiple traits
    bins <- paste0("groupingBIN",1:4)
    # rg
    rg_matrix_birthyear <- all_data %>%
                              filter(Trait1 == PHENOTYPE & Trait2 %in% bins) %>%
                              distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                              select(Trait1, Trait2, rg) %>%
                              pivot_wider(names_from = Trait2, values_from = rg) %>%
                              column_to_rownames("Trait1")
    # rg_se
    rg_se_matrix_birthyear <- all_data %>%
                              filter(Trait1 == PHENOTYPE & Trait2 %in% bins) %>%
                              distinct(Trait1, Trait2, .keep_all = TRUE) %>%
                              select(Trait1, Trait2, se) %>%
                              pivot_wider(names_from = Trait2, values_from = se) %>%
                              column_to_rownames("Trait1")

    # ---------------- Output ----------------
    # Print matrix
    #print(rg_matrix)

    # Save results in a list: all_matrices
    rg_matrix_df <- as.data.frame(rg_matrix_birthyear)
    rg_se_matrix_df <- as.data.frame(rg_se_matrix_birthyear)
    #rg_matrix_df$Phenotype <- PHENOTYPE # keep phenotype label
    all_matrices[[PHENOTYPE]] <- rg_matrix_df
    all_matrices_se[[PHENOTYPE]] <- rg_se_matrix_df

    } #for: PHENOTYPE

    # --------------------------------
    # Combine phenotypes for this GWAS
    # --------------------------------
    # rg
    cat("\nFinal Genetic Correaltion Matrix:\n")
    mtx_name <- paste("GC_", sample,"_", GWAS)
    final_df <- dplyr::bind_rows(all_matrices)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name, final_df)
    print(get(mtx_name))
    # rg_se
    cat("\nFinal SE of Genetic Correaltion Matrix:\n")
    mtx_name_se <- paste("GC_se_", sample,"_", GWAS)
    final_df_se <- dplyr::bind_rows(all_matrices_se)
    #colnames(final_df)[colnames(final_df) == "Smoking"] <- sample
    assign(mtx_name_se, final_df_se)
    print(get(mtx_name_se))
    
    # Store directly
    all_results[[sample]][[GWAS]] <- get(mtx_name)
    all_results_se[[sample]][[GWAS]] <- get(mtx_name_se)
           
} #for: GWAS

# ----------------- rg -------------------- #
cat("\n# ----------- [rg] ----------- #\n")
data_list <- copy(all_results)

# Okbay et al. (2022): EA
rows_to_select <- c("EA")
Full_data <- copy(data_list[["Full"]][["Okbay_et_al_2022"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Okbay_et_al_2022"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Okbay_et_al_2022"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_okbay_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Okbay et al. (2022)]")
combined_okbay <- combined_okbay_all %>%
                    slice(match(rows_to_select, rownames(combined_okbay_all)))
combined_okbay

# Becker et al. (2021): RISK, FAMSAT, FINSAT, FRIENDSAT, WORKSAT
rows_to_select <- c("RISK", "FAMSAT", "FINSAT", "FRIENDSAT", "WORKSAT")
Full_data <- copy(data_list[["Full"]][["Becker_et_al_2021"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Becker_et_al_2021"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Becker_et_al_2021"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_becker_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Bekcer et al. (2021)]")
combined_becker <- combined_becker_all %>%
                    slice(match(rows_to_select, rownames(combined_becker_all)))
combined_becker

# NEALE v2: TDI, RELIGATT
rows_to_select <- c("TDI", "RELIGATT")
Full_data <- copy(data_list[["Full"]][["NEALE_v2"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["NEALE_v2"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["NEALE_v2"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_neale_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[NEALE_v2]")
combined_neale <- combined_neale_all %>%
                        slice(match(rows_to_select, rownames(combined_neale_all)))
combined_neale

# combine data
combined_data <- rbind(combined_okbay, combined_becker, combined_neale)
combined_data
#combined_data$Phenotype <- rownames(combined_data)
column_order <- c("Full", "Female", "Male", "Gen1", "Gen2", "Gen3", "Gen4")
colnames(combined_data) <- column_order

causes_data <-copy(combined_data)

# ----------------- rg_se -------------------- #
cat("\n# ----------- [rg_se] ----------- #\n")
data_list <- copy(all_results_se)

# Okbay et al. (2022): EA
rows_to_select <- c("EA")
Full_data <- copy(data_list[["Full"]][["Okbay_et_al_2022"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Okbay_et_al_2022"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Okbay_et_al_2022"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_okbay_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Okbay et al. (2022)]")
combined_okbay <- combined_okbay_all %>%
                    slice(match(rows_to_select, rownames(combined_okbay_all)))
combined_okbay

# Becker et al. (2021): RISK, FAMSAT, FINSAT, FRIENDSAT, WORKSAT
rows_to_select <- c("RISK", "FAMSAT", "FINSAT", "FRIENDSAT", "WORKSAT")
Full_data <- copy(data_list[["Full"]][["Becker_et_al_2021"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Becker_et_al_2021"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Becker_et_al_2021"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_becker_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Bekcer et al. (2021)]")
combined_becker <- combined_becker_all %>%
                    slice(match(rows_to_select, rownames(combined_becker_all)))
combined_becker

# NEALE v2: TDI, RELIGATT
rows_to_select <- c("TDI", "RELIGATT")
Full_data <- copy(data_list[["Full"]][["NEALE_v2"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["NEALE_v2"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["NEALE_v2"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_neale_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[NEALE_v2]")
combined_neale <- combined_neale_all %>%
                        slice(match(rows_to_select, rownames(combined_neale_all)))
combined_neale

# combine data
combined_data <- rbind(combined_okbay, combined_becker, combined_neale)
combined_data
#combined_data$Phenotype <- rownames(combined_data)
column_order <- c("Full", "Female", "Male", "Gen1", "Gen2", "Gen3", "Gen4")
colnames(combined_data) <- column_order

causes_data_se <-copy(combined_data)

# ----------------- rg -------------------- #
cat("\n# ----------- [rg] ----------- #\n")
data_list <- copy(all_results)

# Liu et al. (2019): SI, CPD, AGESI, SC, DPW
rows_to_select <- c("SI", "CPD", "AGESI", "SC", "DPW")
Full_data <- copy(data_list[["Full"]][["Liu_et_al_2019"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Liu_et_al_2019"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Liu_et_al_2019"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_liu_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Liu et al. (2019)]")
combined_liu <- combined_liu_all %>%
                    slice(match(rows_to_select, rownames(combined_liu_all)))
combined_liu

# Becker et al. (2021): CANNABIS
rows_to_select <- c("CANNABIS")
Full_data <- copy(data_list[["Full"]][["Becker_et_al_2021"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Becker_et_al_2021"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Becker_et_al_2021"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_becker_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Becker et al. (2021)]")
combined_becker <- combined_becker_all %>%
                        slice(match(rows_to_select, rownames(combined_becker_all)))
combined_becker

# combine data
combined_data <- rbind(combined_liu, combined_becker)
combined_data
#combined_data$Phenotype <- rownames(combined_data)
column_order <- c("Full", "Female", "Male", "Gen1", "Gen2", "Gen3", "Gen4")
colnames(combined_data) <- column_order

correlates_data <-copy(combined_data)

# ----------------- rg_se -------------------- #
cat("\n# ----------- [rg_se] ----------- #\n")
data_list <- copy(all_results_se)

# Liu et al. (2019): SI, CPD, AGESI, SC, DPW
rows_to_select <- c("SI", "CPD", "AGESI", "SC", "DPW")
Full_data <- copy(data_list[["Full"]][["Liu_et_al_2019"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Liu_et_al_2019"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Liu_et_al_2019"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_liu_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Liu et al. (2019)]")
combined_liu <- combined_liu_all %>%
                    slice(match(rows_to_select, rownames(combined_liu_all)))
combined_liu

# Becker et al. (2021): CANNABIS
rows_to_select <- c("CANNABIS")
Full_data <- copy(data_list[["Full"]][["Becker_et_al_2021"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Becker_et_al_2021"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Becker_et_al_2021"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_becker_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Becker et al. (2021)]")
combined_becker <- combined_becker_all %>%
                        slice(match(rows_to_select, rownames(combined_becker_all)))
combined_becker

# combine data
combined_data <- rbind(combined_liu, combined_becker)
combined_data
#combined_data$Phenotype <- rownames(combined_data)
column_order <- c("Full", "Female", "Male", "Gen1", "Gen2", "Gen3", "Gen4")
colnames(combined_data) <- column_order

correlates_data_se <-copy(combined_data)

# ----------------- rg -------------------- #
cat("\n# ----------- [rg] ----------- #\n")
data_list <- copy(all_results)

# Becker et al. (2021): BMI
rows_to_select <- c("BMI")
Full_data <- copy(data_list[["Full"]][["Becker_et_al_2021"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Becker_et_al_2021"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Becker_et_al_2021"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_becker_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Becker et al. (2021)]")
combined_becker <- combined_becker_all %>%
                        slice(match(rows_to_select, rownames(combined_becker_all)))
combined_becker

# NEALE v2: ASTHMA, COPD, HAYFEVER, HDL, DIABII, STROKE, HAMI, LUNG, HYPTENS, LONELY
rows_to_select <- c("ASTHMA", "COPD", "HAYFEVER", "HDL", "DIABII", "STROKE", 
                    "HAMI", "LUNG", "HYPTENS", "LONELY")
Full_data <- copy(data_list[["Full"]][["NEALE_v2"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["NEALE_v2"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["NEALE_v2"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_neale_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[NEALE_v2]")
combined_neale <- combined_neale_all %>%
                        slice(match(rows_to_select, rownames(combined_neale_all)))
combined_neale

# Gupta et al. (2024): NEURO, AGREE, CONSC, OPEN, EXTRA
rows_to_select <- c("NEURO", "AGREE", "CONSC", "OPEN", "EXTRA")
Full_data <- copy(data_list[["Full"]][["Gupta_et_al_2024"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Gupta_et_al_2024"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Gupta_et_al_2024"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_gupta_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Gupta et al. (2024)]")
combined_gupta <- combined_gupta_all %>%
                        slice(match(rows_to_select, rownames(combined_gupta_all)))
combined_gupta

# PGC: SCHIZ, ADHD, MDD, ASD, BIP
rows_to_select <- c("SCHIZ", "ADHD", "MDD", "ASD", "BIP")
Full_data <- copy(data_list[["Full"]][["PGC"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["PGC"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["PGC"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_pgc_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[PGC]")
combined_pgc <- combined_pgc_all %>%
                        slice(match(rows_to_select, rownames(combined_pgc_all)))
combined_pgc


# combine data
combined_data <- rbind(combined_becker, combined_neale, combined_gupta, combined_pgc)
combined_data
#combined_data$Phenotype <- rownames(combined_data)
column_order <- c("Full", "Female", "Male", "Gen1", "Gen2", "Gen3", "Gen4")
colnames(combined_data) <- column_order

consequences_data <-copy(combined_data)

# ----------------- rg_se -------------------- #
cat("\n# ----------- [rg_se] ----------- #\n")
data_list <- copy(all_results_se)

# Becker et al. (2021): BMI
rows_to_select <- c("BMI")
Full_data <- copy(data_list[["Full"]][["Becker_et_al_2021"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Becker_et_al_2021"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Becker_et_al_2021"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_becker_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Becker et al. (2021)]")
combined_becker <- combined_becker_all %>%
                        slice(match(rows_to_select, rownames(combined_becker_all)))
combined_becker

# NEALE v2: ASTHMA, COPD, HAYFEVER, HDL, DIABII, STROKE, HAMI, LUNG, HYPTENS, LONELY
rows_to_select <- c("ASTHMA", "COPD", "HAYFEVER", "HDL", "DIABII", "STROKE", 
                    "HAMI", "LUNG", "HYPTENS", "LONELY")
Full_data <- copy(data_list[["Full"]][["NEALE_v2"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["NEALE_v2"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["NEALE_v2"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_neale_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[NEALE_v2]")
combined_neale <- combined_neale_all %>%
                        slice(match(rows_to_select, rownames(combined_neale_all)))
combined_neale

# Gupta et al. (2024): NEURO, AGREE, CONSC, OPEN, EXTRA
rows_to_select <- c("NEURO", "AGREE", "CONSC", "OPEN", "EXTRA")
Full_data <- copy(data_list[["Full"]][["Gupta_et_al_2024"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["Gupta_et_al_2024"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["Gupta_et_al_2024"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_gupta_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[Gupta et al. (2024)]")
combined_gupta <- combined_gupta_all %>%
                        slice(match(rows_to_select, rownames(combined_gupta_all)))
combined_gupta

# PGC: SCHIZ, ADHD, MDD, ASD, BIP
rows_to_select <- c("SCHIZ", "ADHD", "MDD", "ASD", "BIP")
Full_data <- copy(data_list[["Full"]][["PGC"]])
Full_data[] <- lapply(Full_data, as.numeric) 
gender_data <- copy(data_list[["gender"]][["PGC"]])
gender_data[] <- lapply(gender_data, as.numeric) 
generation_data <- copy(data_list[["generation"]][["PGC"]])
generation_data[] <- lapply(generation_data, as.numeric) 
# generate combined data: column bind
combined_pgc_all <- cbind(Full_data,
                           gender_data,
                           generation_data
                          )
print("[PGC]")
combined_pgc <- combined_pgc_all %>%
                        slice(match(rows_to_select, rownames(combined_pgc_all)))
combined_pgc


# combine data
combined_data <- rbind(combined_becker, combined_neale, combined_gupta, combined_pgc)
combined_data
#combined_data$Phenotype <- rownames(combined_data)
column_order <- c("Full", "Female", "Male", "Gen1", "Gen2", "Gen3", "Gen4")
colnames(combined_data) <- column_order

consequences_data_se <-copy(combined_data)

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

# replace rownames
## [causes_data]
rownames(causes_data) <- ifelse(
                          rownames(causes_data) %in% names(pheno_vec),
                          pheno_vec[rownames(causes_data)],
                          rownames(causes_data)
                        )
## [causes_data_se]
rownames(causes_data_se) <- ifelse(
                              rownames(causes_data_se) %in% names(pheno_vec),
                              pheno_vec[rownames(causes_data_se)],
                              rownames(causes_data_se)
                            )

## [correlates_data]
rownames(correlates_data) <- ifelse(
                              rownames(correlates_data) %in% names(pheno_vec),
                              pheno_vec[rownames(correlates_data)],
                              rownames(correlates_data)
                            )
## [correlates_data_se]
rownames(correlates_data_se) <- ifelse(
                              rownames(correlates_data_se) %in% names(pheno_vec),
                              pheno_vec[rownames(correlates_data_se)],
                              rownames(correlates_data_se)
                            )

## [consequences_data]
rownames(consequences_data) <- ifelse(
                              rownames(consequences_data) %in% names(pheno_vec),
                              pheno_vec[rownames(consequences_data)],
                              rownames(consequences_data)
                            )
## [consequences_data_se]
rownames(consequences_data_se) <- ifelse(
                              rownames(consequences_data_se) %in% names(pheno_vec),
                              pheno_vec[rownames(consequences_data_se)],
                              rownames(consequences_data_se)
                            )


#combine data
cat("\n[data]\n")
data <- rbind(causes_data, correlates_data, consequences_data)
data
cat("\n[data_se]\n")
data_se <- rbind(causes_data_se, correlates_data_se, consequences_data_se)
data_se


# -------------------------
# 1) Block labels
# -------------------------
causes_names <- rownames(causes_data)
correlates_names <- rownames(correlates_data)
consequences_names <- rownames(consequences_data)

# Preserve original order
pheno_order <- c(rownames(causes_data), rownames(correlates_data),rownames(consequences_data))
pheno_order <- rev(pheno_order) # top-to-bottom order

data$Block <- NA
data$Block[rownames(data) %in% causes_names] <- "Sociodemographics"
data$Block[rownames(data) %in% correlates_names] <- "Substance Use"
data$Block[rownames(data) %in% consequences_names] <- "Psychiatric and Medical Conditions"

data$Phenotype <- rownames(data)
data
# rename the 7 columns (adjust if your names differ)
column_names <- c("Full","Female","Male","Gen1","Gen2","Gen3","Gen4")
colnames(data)[1:7] <- column_names

# make facet order explicit
data$Block <- factor(data$Block, levels = c("Substance Use","Sociodemographics","Psychiatric and Medical Conditions"))


# -------------------------
# 2) Long format
# -------------------------
# rg
df_long <- data %>%
  pivot_longer(
    cols = all_of(column_names),
    names_to = "Group",
    values_to = "rg"
  ) %>%
  mutate(
    rg = as.numeric(rg),
    Group = factor(Group, levels = column_names)
  ) #%>%
  #group_by(Block) %>%
  #mutate(
    # IMPORTANT: facet-specific phenotype factor (prevents huge blank space)
  #  Phenotype = factor(Phenotype, levels = rev(unique(Phenotype)))
  #) %>%
  #ungroup()

# rg_se
df_se_long <- data_se %>%
  mutate(
    Block = data$Block,
    Phenotype = rownames(data_se)
  ) %>%
  pivot_longer(
    cols = all_of(column_names),
    names_to = "Group",
    values_to = "SE"
  ) %>%
  mutate(
    Group = factor(Group, levels = column_names),
    SE = as.numeric(SE)
  )

## Merge [rg] and [rg_se]
df_test <- df_long %>%
  left_join(df_se_long, by = c("Phenotype","Group","Block"))

# Apply fixed factor levels after merging
df_test$Phenotype <- factor(df_test$Phenotype, levels = pheno_order)

# -------------------------
# 3) Test for Significance ([rg=0])
#    with Group specific bonferroni correction.
# -------------------------
cat("\nNumber of Tests: \n")
n_pheno <- length(unique(df_test$Phenotype))
n_pheno

df_test <- df_test %>%
  group_by(Group) %>%
  mutate(
    num_tests = case_when (
        (Group == "Full") ~ (1*n_pheno),
        (Group %in% c("Female", "Male")) ~ (2*n_pheno),
        (Group %in% c("Gen1", "Gen2", "Gen3", "Gen4")) ~ (4*n_pheno)
    ),
    z = rg / SE,
    p = 2 * pnorm(-abs(z)),
    stars = case_when(
      p < (0.05 / num_tests) ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

cat("\n[df_test]\n")
#df_test
cat("\n[table(df_test$Group)]\n")
table(df_test$Group)
cat("\n[Group is Full]\n")
df_test[df_test$Group=="Full",]


# -------------------------
# 5) Plot
# -------------------------
gc_traits <- ggplot(df_test, aes(Group, Phenotype, fill = rg)) +
              geom_tile(color = "white", linewidth = 0.7) +
              geom_text(aes(label = paste0(sprintf("%.2f", rg), stars)), size = 3.2) +
              geom_vline(xintercept = c(1.5, 3.5), color = "white", linewidth = 2.5) +
              scale_fill_gradient2(
              low = "#A020F0",   # bright purple
              mid = "white",
              high = "#FF0000",  # bright red
                midpoint = 0,
                limits = c(-1, 1),
                name = "Genetic Correlation",
                na.value = "white"
              ) +
              facet_grid(Block ~ ., scales = "free_y", space = "free_y") +
                scale_x_discrete(position = "top",
                  labels = c("Full" = "Full\nSample",
                             "Female" = "\nFemale",
                             "Male" = "\nMale",
                             "Gen1" = "Silent\nGeneration",
                             "Gen2" = "Baby\nBoomers",
                             "Gen3" = "Generation\nX",
                             "Gen4" = "\nMillennials")) +
              scale_y_discrete(position = "right") +
              theme_minimal(base_size = 12) +
              theme(
                panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text.x.top = element_text(size = 8, lineheight = 0.85),
                axis.text.y.right = element_text(hjust = 0),
                strip.background.y = element_rect(fill = "grey90", color = "black"),
                strip.text.y = element_text(angle = 90, face = "bold", size = 7),
                legend.position = "none"
              )
gc_traits

# Save it
ggsave(paste0(working_dir, "/output/Figures/rg_all_traits.png"), plot=gc_traits, width = 10, height = 6, dpi = 500)

df_test %>%
 filter(Phenotype == "Lung Cancer")

# -----------------------------------------
# Pairwise tests of equality (rg_A==rg_B)
# -----------------------------------------

# ---------------------------------------------------------------
# 1) Define All Pairwise Comparisons, sample size, and pooled SE
# ---------------------------------------------------------------
pairwise_list <- list(
  "Female vs Male" = c("Female","Male"),
  "Gen1 vs Gen2"   = c("Gen1","Gen2"),
  "Gen1 vs Gen3"   = c("Gen1","Gen3"),
  "Gen1 vs Gen4"   = c("Gen1","Gen4"),
  "Gen2 vs Gen3"   = c("Gen2","Gen3"),
  "Gen2 vs Gen4"   = c("Gen2","Gen4"),
  "Gen3 vs Gen4"   = c("Gen3","Gen4")
)

# Number of sample sizes
sample_sizes <- c(
    Female = 2440320,
    Male = 2016119,
    Gen1 = 23435 + 50914 + 99204 + 192470,
    Gen2 = 328509 + 415050 + 459411 + 474935,
    Gen3 = 461498 + 450003 + 419416 + 495137,
    Gen4 = 541787 + 511499 + 356082 + 190406
) #sample_sizes

# define function for [pooled se]
pooled_se <- function(n1, se1, n2, se2){
  pooled <- sqrt((((n1 - 1)*se1^2) + ((n2 - 1)*se2^2)) / (n1 + n2 - 2))
  pooled
}

# -----------------------------------------
# 2) Compute Pairwise Difference Tests
# -----------------------------------------
pairwise_results <- lapply(names(pairwise_list), function(comp){

  g1 <- pairwise_list[[comp]][1]
  g2 <- pairwise_list[[comp]][2]

  df_test %>%
    filter(Group %in% c(g1, g2)) %>%
    select(Phenotype, Block, Group, rg, SE) %>%
    pivot_wider(names_from = Group,
                values_from = c(rg, SE)) %>%
    mutate(
      # difference
      diff = get(paste0("rg_", g1)) -
             get(paste0("rg_", g2)),      
      se_diff = pooled_se(
                  sample_sizes[g1],
                  get(paste0("SE_", g1)),
                  sample_sizes[g2],
                  get(paste0("SE_", g2))
                  ),      
      z = diff / se_diff,
      p = 2 * pnorm(-abs(z)),
      
      Comparison = comp
    ) %>%
    select(Block, Phenotype, Comparison, diff, z, p)
})

pairwise_results <- bind_rows(pairwise_results)
#pairwise_results

# -----------------------------------------
# 3) Add significance
# -----------------------------------------
n_pheno <- length(unique(pairwise_results$Phenotype))

pairwise_results <- pairwise_results %>%
  mutate(
    family_size = case_when(
      (Comparison == "Female vs Male") ~ n_pheno,
      TRUE ~ 6 * n_pheno
    ),
    p_bonf = p * family_size,
    p_bonf = pmin(p_bonf, 1),

    symbol = case_when(
      p_bonf < 0.05 ~ "⏺⏺",   # Bonferroni significant
      p < 0.05 ~ "⏺",         # Nominal significant
      TRUE ~ "⭘"
    )
  )

cat("[pairwise_results]: ")
#pairwise_results

# -----------------------------------------
# 4) Pivot to Wide Table
# -----------------------------------------
table_output <- pairwise_results %>%
  select(Block, Phenotype, Comparison, symbol) %>%
  pivot_wider(names_from = Comparison,
              values_from = symbol)

cat("[table_output]: ")
table_output


# -----------------------------------------
# 5) Export as Excel Table
# -----------------------------------------
cat("\nExport [table_output]: \n")
#gt_table <- gt(table_output)
#gtsave(gt_table, paste0(working_dir, "/output/Tables/pairwise_results_table.html"))
write_xlsx(table_output, paste0(working_dir, "/output/Tables/pairwise_results_table.xlsx"))
cat("\n* --------------------------------------------- *\n")

# -------------------------
# 1) Block labels
# -------------------------
causes_names <- rownames(causes_data)
correlates_names <- rownames(correlates_data)
consequences_names <- rownames(consequences_data)

# Preserve original order
pheno_order <- c(rownames(causes_data), rownames(correlates_data),rownames(consequences_data))
pheno_order <- rev(pheno_order) # top-to-bottom order

data$Block <- NA
data$Block[rownames(data) %in% causes_names] <- "Sociodemographics"
data$Block[rownames(data) %in% correlates_names] <- "Substance Use"
data$Block[rownames(data) %in% consequences_names] <- "Psychiatric and Medical Conditions"

data$Phenotype <- rownames(data)

# rename the 7 columns (adjust if your names differ)
column_names <- c("Full","Female","Male","Gen1","Gen2","Gen3","Gen4")
colnames(data)[1:7] <- column_names

# make facet order explicit
data$Block <- factor(data$Block, levels = c("Substance Use","Sociodemographics","Psychiatric and Medical Conditions"))


# -------------------------
# 2) Long format
# -------------------------
# rg
df_long <- data %>%
  pivot_longer(
    cols = all_of(column_names),
    names_to = "Group",
    values_to = "rg"
  ) %>%
  mutate(
    rg = as.numeric(rg),
    Group = factor(Group, levels = column_names)
  ) #%>%
  #group_by(Block) %>%
  #mutate(
    # IMPORTANT: facet-specific phenotype factor (prevents huge blank space)
  #  Phenotype = factor(Phenotype, levels = rev(unique(Phenotype)))
  #) %>%
  #ungroup()

# rg_se
df_se_long <- data_se %>%
  mutate(
    Block = data$Block,
    Phenotype = rownames(data_se)
  ) %>%
  pivot_longer(
    cols = all_of(column_names),
    names_to = "Group",
    values_to = "SE"
  ) %>%
  mutate(
    Group = factor(Group, levels = column_names),
    SE = as.numeric(SE)
  )

## Merge [rg] and [rg_se]
df_test <- df_long %>%
  left_join(df_se_long, by = c("Phenotype","Group","Block"))

# Apply fixed factor levels after merging
df_test$Phenotype <- factor(df_test$Phenotype, levels = pheno_order)
#df_test

# -------------------------
# 3) Test for Significance ([rg=0])
#    with Group specific bonferroni correction.
# -------------------------
cat("\nNumber of Tests: \n")
n_pheno <- length(unique(df_test$Phenotype))
n_pheno

df_test <- df_test %>%
  group_by(Group) %>%
  mutate(
    num_tests = case_when (
        (Group == "Full") ~ (1*n_pheno),
        (Group %in% c("Female", "Male")) ~ (2*n_pheno),
        (Group %in% c("Gen1", "Gen2", "Gen3", "Gen4")) ~ (4*n_pheno)
    ),
    z = rg / SE,
    p = 2 * pnorm(-abs(z)),
    stars = case_when(
      p < (0.05 / num_tests) ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

cat("\n[df_test]]\n")
#df_test


P <- length(unique(df_test$Phenotype))
# -----------------------------------------
# 4) Female vs Male difference test
# -----------------------------------------
df_fm <- df_test %>%
  filter(Group %in% c("Female","Male")) %>%
  select(Phenotype, Block, Group, rg, SE) %>%
  pivot_wider(names_from = Group,
              values_from = c(rg, SE)) %>%
  mutate(
    # difference
    diff = (rg_Male - rg_Female),
    se_diff  = pooled_se(sample_sizes["Female"], 
                         get(paste0("SE_Female")),
                         sample_sizes["Male"],
                         get(paste0("SE_Male"))),
    z        = diff / se_diff,
    p        = 2 * pnorm(-abs(z)),

    # Bonferroni = 1 × P
    p_bonf   = p * P,
    p_bonf   = pmin(p_bonf, 1),

    symbol = case_when(
      p_bonf < 0.05 ~ "**",
      p < 0.05      ~ "*",
      TRUE          ~ ""
    ),

    Group = "F vs M"
  ) %>%
  select(Phenotype, Block, Group, symbol)
cat("\n[df_fm]\n")
df_fm

# --------------------------------------------------------------------------------------
# 5) (Gen1=Gen2), (Gen1=Gen3), (Gen1=Gen4), (Gen2=Gen3), (Gen2=Gen4), (Gen3=Gen4) tests
# --------------------------------------------------------------------------------------
#Extract Gen1–Gen4
df_gen <- df_test %>%
  filter(Group %in% c("Gen1","Gen2","Gen3","Gen4")) %>%
  select(Phenotype, Block, Group, rg, SE) %>%
  pivot_wider(names_from = Group,
              values_from = c(rg, SE))

#Define all 6 pairwise comparisons
gen_pairs <- list(
  "Gen1 vs Gen2" = c("Gen1","Gen2"),
  "Gen1 vs Gen3" = c("Gen1","Gen3"),
  "Gen1 vs Gen4" = c("Gen1","Gen4"),
  "Gen2 vs Gen3" = c("Gen2","Gen3"),
  "Gen2 vs Gen4" = c("Gen2","Gen4"),
  "Gen3 vs Gen4" = c("Gen3","Gen4")
) #gen_pairs

#Compute Wald tests for each pair
df_gen_test <- lapply(names(gen_pairs), function(comp){

  g1 <- gen_pairs[[comp]][1]
  g2 <- gen_pairs[[comp]][2]

  df_gen %>%
    mutate(
      diff = get(paste0("rg_", g1)) -
             get(paste0("rg_", g2)),

      se_diff = pooled_se(sample_sizes[g1],
                          get(paste0("SE_", g1)),
                          sample_sizes[g2],
                          get(paste0("SE_", g2))
                         ),
      z = diff / se_diff,
      p = 2 * pnorm(-abs(z)),

      Comparison = comp
    ) %>%
    select(Block, Phenotype, Comparison, p)
})

df_gen_test <- bind_rows(df_gen_test)


# Bonferroni correction (Num of test = 6 × P)
n_pheno <- length(unique(df_gen_test$Phenotype))

df_gen_test <- df_gen_test %>%
  mutate(
    p_bonf = p * (6 * n_pheno),
    p_bonf = pmin(p_bonf, 1),

    symbol = case_when(
      p_bonf < 0.05 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ), #case_when
    Group = "Gen Joint"
  ) #%>%
  #select(Phenotype, Block, Group, symbol)
cat("\n[df_gen_test]\n")
#df_gen_test

# Make wide table of pairwise results
df_gen_wide <- df_gen_test %>%
  select(Block, Phenotype, Comparison, symbol, p, p_bonf) %>%
  pivot_wider(
    names_from = Comparison,
    values_from = c(symbol, p, p_bonf)
  )

cat("\n[df_gen_wide]\n")
#df_gen_wide

# Create Gen Joint column
p_cols <- grep("^p_", names(df_gen_wide), value = TRUE)
pbonf_cols <- grep("^p_bonf_", names(df_gen_wide), value = TRUE)

df_gen_wide <- df_gen_wide %>%
  rowwise() %>%
  mutate(
    any_bonf = any(c_across(all_of(pbonf_cols)) < 0.05, na.rm = TRUE),
    any_nominal = any(c_across(all_of(p_cols)) < 0.05, na.rm = TRUE),

    Gen_Joint = case_when(
      any_bonf ~ "**",
      any_nominal ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

cat("\n[df_gen_wide with Gen_Joint]\n")
#df_gen_wide

# Select necessary columns
df_gen_final <- df_gen_wide %>%
                  select(Block, Phenotype, Gen_Joint) %>%
                  rename(symbol = Gen_Joint) %>%
                  mutate(Group = "Gen Joint")
cat("\n[df_gen_final]\n")
df_gen_final


#Add new columns
df_plot2 <- df_test %>%
  select(Phenotype, Block, Group, rg, stars) %>%
  bind_rows(
    df_fm %>% mutate(rg = NA),
    df_gen_final %>% mutate(rg = NA)
  )

df_plot2$Group <- factor(
  df_plot2$Group,
  levels = c("Full","Female","Male","F vs M",
             "Gen1","Gen2","Gen3","Gen4","Gen Joint")
)

cat("\n[df_plot2]\n")
df_plot2

# create dummy rows for spacer column
df_spacer <- df_plot2 %>%
  distinct(Phenotype, Block) %>%
  mutate(
    Group = "Spacer",
    rg = NA,
    stars = "",
    symbol = ""
  )

# combine with original data
df_plot2_new <- bind_rows(df_plot2, df_spacer)

# set desired order
df_plot2_new$Group <- factor(
  df_plot2_new$Group,
  levels = c(
    "Full", "Spacer", "Female", "Male", "F vs M",
    "Gen1", "Gen2", "Gen3", "Gen4", "Gen Joint"
  )
)

# -------------------------
# 6) Plot
# -------------------------
gc_traits_new <- ggplot(df_plot2_new, aes(Group, Phenotype, fill = rg)) +
              geom_tile(color = "white", linewidth = 0.7) +
              geom_text(aes(label = ifelse(Group %in% c("F vs M", "Gen Joint", "Spacer"), symbol, 
                                           paste0(sprintf("%.2f", rg), stars))), size = 3.2) +
              scale_fill_gradient2(
              low = "#A020F0",   # bright purple
              mid = "white",
              high = "#FF0000",  # bright red
                midpoint = 0,
                limits = c(-1, 1),
                name = "Genetic Correlation",
                na.value = "grey90"
              ) +
              facet_grid(Block ~ ., scales = "free_y", space = "free_y") +
              scale_x_discrete(position = "top",
                               labels = c("Full" = "Full\nSample",
                                          "Spacer" = "",
                                          "Female" = "\nFemale",
                                          "Male" = "\nMale",
                                          "F vs M" = "\n\u03b4sex",
                                          "Gen1" = "Silent\nGeneration",
                                          "Gen2" = "Baby\nBoomers",
                                          "Gen3" = "Generation\nX",
                                          "Gen4" = "\nMillennials",
                                          "Gen Joint" = "\n\u03b4gen")) +        
              scale_y_discrete(position = "right") +
              theme_minimal(base_size = 12) +
              theme(
                panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text.x.top = element_text(size = 8,lineheight = 0.85),
                axis.text.y.right = element_text(hjust = 0),
                strip.background.y = element_rect(fill = "grey90", color = "black"),
                strip.text.y = element_text(angle = 90, face = "bold", size = 7),
                legend.position = "none"
              )
gc_traits_new

#Save
ggsave(paste0(working_dir, "/output/Figures/rg_all_traits_new.png"), plot=gc_traits_new, width = 10, height = 6, dpi = 500)

#Current directory
cat("Current directory\n")
setwd(root_dir)
cat(getwd())
cat("\n-------------------------------------------------------------------------------\n")

# set some directories
heritability_SI_dir=paste0(root_dir,"/Heritability/raw")

# Function to extract heritability summary from each log file
extract_h2_summary <- function(file) {

  lines <- readLines(file)
  
# ---- Extract sample name from file path ----
  if (grepl("ever_tobacco_user_[a-zA-Z]+_h2", file)) {
    sample <- sub(".*ever_tobacco_user_([a-zA-Z]+)_h2.*", "\\1", file)
    
  } else if (grepl("birthyear_bin[0-9]+", file)) {
    sample <- sub(".*birthyear_bin([0-9]+).*", "bin\\1", file)
    
  } else if (grepl("groupingBIN[0-9]+", file)) {
    sample <- sub(".*groupingBIN([0-9]+).*", "groupingBIN\\1", file)
    
  } else {
    sample <- "full"
  }
    
# ---- Extract LDSC output lines ----
  h2_line         <- grep("Total Observed scale h2:", lines, value = TRUE)
  intercept_line  <- grep("^Intercept:", lines, value = TRUE)
  lambda_line     <- grep("^Lambda GC:", lines, value = TRUE)
  mean_chi2_line  <- grep("^Mean Chi\\^2:", lines, value = TRUE)
  ratio_line      <- grep("^Ratio:", lines, value = TRUE)

# ---- Parse numeric values safely ----
  h2          <- if (length(h2_line)) as.numeric(sub(".*h2: ([0-9.]+).*", "\\1", h2_line)) else NA
  h2_se       <- if (length(h2_line)) as.numeric(sub(".*\\((.*)\\).*", "\\1", h2_line)) else NA
  intercept   <- if (length(intercept_line)) as.numeric(sub(".*Intercept: ([0-9.]+).*", "\\1", intercept_line)) else NA
  intercept_se<- if (length(intercept_line)) as.numeric(sub(".*\\((.*)\\).*", "\\1", intercept_line)) else NA
  lambda_gc   <- if (length(lambda_line)) as.numeric(sub(".*Lambda GC: ([0-9.]+).*", "\\1", lambda_line)) else NA
  mean_chi2   <- if (length(mean_chi2_line)) as.numeric(sub(".*Mean Chi\\^2: ([0-9.]+).*", "\\1", mean_chi2_line)) else NA
  ratio       <- if (length(ratio_line)) as.numeric(sub(".*Ratio: ([0-9.]+).*", "\\1", ratio_line)) else NA

# ---- Return one-row data frame ----
  data.frame(
    Sample = sample,
    h2 = h2,
    h2_se = h2_se,
    Intercept = intercept,
    Intercept_se = intercept_se,
    Lambda_GC = lambda_gc,
    Mean_Chi2 = mean_chi2,
    Ratio = ratio
  )
}


# ---------------------------- FULL sample -------------------------#

# Target subfolders
subdirs <- c("ever_tobacco_user")

# Construct full log file paths
log_files <- file.path(heritability_SI_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# Rename Sample
heritability_df$Sample[heritability_df$Sample == "full"] <- "Full"

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# assign
heritability_full <- copy(heritability_df)

# View results
print(heritability_full)

# ---------------------------- FULL and GENDER sample -------------------------#

# Target subfolders
subdirs <- c("ever_tobacco_user_female", "ever_tobacco_user_male")

# Construct full log file paths
log_files <- file.path(heritability_SI_dir, paste0(subdirs, "_h2_raw.log"))


# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# Rename Sample
heritability_df$Sample[heritability_df$Sample == "female"] <- "Female"
heritability_df$Sample[heritability_df$Sample == "male"] <- "Male"

# assign
heritability_gender <- copy(heritability_df)

# View results
print(heritability_gender)

# ------------------------------- REGION sample -------------------------------#

# Target subfolders
subdirs <- paste0("ever_tobacco_user_", c("MW","NE","SE","SW","W"))

# Construct full log file paths
log_files <- file.path(heritability_SI_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# assign
heritability_region <- copy(heritability_df)

# View results
print(heritability_region)

# ----------------------------- BIRTHYEAR sample ------------------------------#

## create a birtyear ranges
 start_byears <- seq(1924, 1999, by=5)
 end_years <- start_byears + 4
 byear_ranges <- paste0(start_byears, "-", end_years)
  print(byear_ranges) 
  typeof(byear_ranges)

# Target subfolders
subdirs <- paste0("ever_tobacco_user_birthyear_bin", 4:19)

# Construct full log file paths
log_files <- file.path(heritability_SI_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# add birthyear range
heritability_df$byear_ranges <- byear_ranges

# assign
heritability_birthyear <- copy(heritability_df)

# View results
print(heritability_birthyear)

# ----------------------------- BIRTHYEAR sample ------------------------------#

## create a birtyear ranges
 start_byears <- seq(1924, 1984, by=20)
 end_years <- start_byears + 19
 byear_ranges <- paste0(start_byears, "-", end_years)
  print(byear_ranges) 
  typeof(byear_ranges)

# Target subfolders
subdirs <- paste0("ever_tobacco_user_birthyear_groupingBIN", 1:4)

# Construct full log file paths
log_files <- file.path(heritability_SI_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# add birthyear range
heritability_df$byear_ranges <- byear_ranges

# assign
heritability_generation <- copy(heritability_df)

# View results
print(heritability_generation)

#Current directory
cat("Current directory\n")
setwd(root_dir)
cat(getwd())
cat("\n-------------------------------------------------------------------------------\n")

# set some directories
gc_SI_dir=paste0(root_dir,"/Genetic_Correlation/raw/SI")

#-------------------------------- Gender --------------------------------------#
# Set working directory to where your LDSC log files are
setwd(file.path(gc_SI_dir, "/ever_tobacco_user_gender"))

# List only relevant log files
log_files <- list.files(pattern = "_rg_raw.log$")

# Function to extract rg summary from one log file
extract_rg_summary <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  if (length(start) == 0) return(NULL)
  data_lines <- lines[(start + 1):length(lines)]
  data_lines <- data_lines[data_lines != ""]
  
  tryCatch({
    df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
    colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                      "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                      "gcov_int", "gcov_int_se")
    
# Extract gender labels ("female", "male") from file paths
    extract_gender <- function(path) {
      if (grepl("female", path)) return("female")
      else if (grepl("male", path)) return("male")
      else return(NA)
    }
    
    df$Gender1 <- sapply(df$p1, extract_gender)
    df$Gender2 <- sapply(df$p2, extract_gender)
    
    return(df)
  }, error = function(e) NULL)
}

# Read and combine all matching logs
all_data <- map_df(log_files, extract_rg_summary)
all_data
           
# Define the region codes of interest
genders <- c("female", "male")

# Filter only desired gender-gender combinations and build matrix
# rg
rg_matrix_gender <- all_data %>%
             filter(Gender1 %in% genders, Gender2 %in% genders) %>%
             distinct(Gender1, Gender2, .keep_all = TRUE) %>%
             select(Gender1, Gender2, rg) %>%
             pivot_wider(names_from = Gender2, values_from = rg) %>%
             column_to_rownames("Gender1") %>%
             select(any_of(genders)) %>%   
             .[intersect(rownames(.), genders), , drop=F] 
# rg_se
rg_se_matrix_gender <- all_data %>%
             filter(Gender1 %in% genders, Gender2 %in% genders) %>%
             distinct(Gender1, Gender2, .keep_all = TRUE) %>%
             select(Gender1, Gender2, se) %>%
             pivot_wider(names_from = Gender2, values_from = se) %>%
             column_to_rownames("Gender1") %>%
             select(any_of(genders)) %>%   
             .[intersect(rownames(.), genders), , drop=F] 


# column and row names
colnames(rg_matrix_gender) <-"Male"
colnames(rg_se_matrix_gender) <-"Male"
rownames(rg_matrix_gender) <- "Female"
rownames(rg_se_matrix_gender) <- "Female"

# Print the matrix
cat("[rg_matrix_gender]\n")
print(rg_matrix_gender)
cat("[rg_se_matrix_gender]\n")
print(rg_se_matrix_gender)

#------------------------------------ Region ---------------------------------------#
# Set working directory to where your LDSC log files are
setwd(file.path(gc_SI_dir, "/ever_tobacco_user_region"))

# List only relevant log files
log_files <- list.files(pattern = "_rg_raw.log$")

# Function to extract rg summary from one log file
extract_rg_summary <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  if (length(start) == 0) return(NULL)
  data_lines <- lines[(start + 1):length(lines)]
  data_lines <- data_lines[data_lines != ""]
  
  tryCatch({
    df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
    colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                      "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                      "gcov_int", "gcov_int_se")
    
    # Extract region names
    extract_region <- function(path) sub(".*/ever_tobacco_user_([A-Z]+)/.*", "\\1", path)
    df$Region1 <- sapply(df$p1, extract_region)
    df$Region2 <- sapply(df$p2, extract_region)
    
    return(df)
  }, error = function(e) NULL)
}

# Read and combine all matching logs
all_data <- map_df(log_files, extract_rg_summary)
all_data
           
# Define the region codes of interest
regions <- c("MW", "NE", "SE", "SW", "W")

# Filter only desired region-region combinations and build matrix
#rg
rg_matrix_region <- all_data %>%
             filter(Region1 %in% regions, Region2 %in% regions) %>%
             distinct(Region1, Region2, .keep_all = TRUE) %>%
             select(Region1, Region2, rg) %>%
             pivot_wider(names_from = Region2, values_from = rg) %>%
             column_to_rownames("Region1") %>%
             select(any_of(regions)) %>%   
             .[intersect(rownames(.), regions), ] 
#rg
rg_se_matrix_region <- all_data %>%
             filter(Region1 %in% regions, Region2 %in% regions) %>%
             distinct(Region1, Region2, .keep_all = TRUE) %>%
             select(Region1, Region2, se) %>%
             pivot_wider(names_from = Region2, values_from = se) %>%
             column_to_rownames("Region1") %>%
             select(any_of(regions)) %>%   
             .[intersect(rownames(.), regions), ] 

# Print the matrix
cat("[rg_matrix_region]\n")
print(rg_matrix_region)
cat("[rg_se_matrix_region]\n")
print(rg_se_matrix_region)

#------------------------------------ Birthyear Bin: bin4 to bin19 ---------------------------------------#

## create a birtyear ranges
 start_byears <- seq(1924, 1999, by=5)
 end_years <- start_byears + 4
 byear_ranges <- paste0(start_byears, "-", end_years)
  print(byear_ranges) 
  typeof(byear_ranges)

## create birthyear_bins
  byear_bins <- paste0("bin",4:19)
  
## combine [birthyear_bins] and  [birthyear_ranges]
  #birthyear <- cbind(byear_bins, byear_ranges)
  birthyear <- list()
  for (i in 1:length(byear_bins)) {
   birthyear[[byear_bins[i]]] <- byear_ranges[i]  # a list
  } #for: i
  
  print(birthyear)
  
# Set working directory to where your LDSC log files are
setwd(file.path(gc_SI_dir, "/ever_tobacco_user_birthyear"))

# List LDSC log files matching the naming pattern
log_files <- list.files(pattern = "9_rg_raw.log$")
log_files

# Function to extract rg summary from one log file
extract_rg_summary <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  if (length(start) == 0) return(NULL)
  data_lines <- lines[(start + 1):length(lines)]
  data_lines <- data_lines[data_lines != ""]
  
  tryCatch({
    df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
    colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                      "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                      "gcov_int", "gcov_int_se")
    
    # Extract bin names in format "bin4", "bin5", ..., "bin19"
    extract_bin <- function(path) sub(".*birthyear_bin(\\d+).*", "bin\\1", path)
    df$Bin1 <- sapply(df$p1, extract_bin)
    df$Bin2 <- sapply(df$p2, extract_bin)
    
    return(df)
  }, error = function(e) NULL)
}


# Extract and combine all data
all_data <- map_df(log_files, extract_rg_summary)

# Define full bin range
bins <- paste0("bin", 4:19)

# Filter, pivot, and clean up matrix
# rg
rg_matrix_birthyear <- all_data %>%
                        filter(Bin1 %in% bins, Bin2 %in% bins) %>%
                        distinct(Bin1, Bin2, .keep_all = TRUE) %>%
                        select(Bin1, Bin2, rg) %>%
                        pivot_wider(names_from = Bin2, values_from = rg) %>%
                        column_to_rownames("Bin1") %>%
                        select(any_of(bins)) %>%
                        .[intersect(rownames(.), bins), ]
# se
rg_se_matrix_birthyear <- all_data %>%
                        filter(Bin1 %in% bins, Bin2 %in% bins) %>%
                        distinct(Bin1, Bin2, .keep_all = TRUE) %>%
                        select(Bin1, Bin2, se) %>%
                        pivot_wider(names_from = Bin2, values_from = se) %>%
                        column_to_rownames("Bin1") %>%
                        select(any_of(bins)) %>%
                        .[intersect(rownames(.), bins), ]


# Print the matrix
cat("\n[rg_matrix_birthyear]\n")
print(rg_matrix_birthyear)
cat("\n[rg_se_matrix_birthyear]\n")
print(rg_se_matrix_birthyear)

# Reorder rows
bins_or <- paste0("bin", 4:18)
rg_matrix_birthyear_ro <- rg_matrix_birthyear[bins_or, ]
rg_se_matrix_birthyear_ro <- rg_se_matrix_birthyear[bins_or, ]

# Print the matrix
cat("\n[rg_matrix_birthyear_ro]\n")
print(rg_matrix_birthyear_ro)
cat("\n[rg_se_matrix_birthyear_ro]\n")
print(rg_se_matrix_birthyear_ro)


#------------------------------------ Birthyear Bin: bin4 to bin19 ---------------------------------------#

## create a birtyear ranges
 start_byears <- seq(1924, 1984, by=20)
 end_years <- start_byears + 19
 byear_ranges <- paste0(start_byears, "-", end_years)
  print(byear_ranges) 
  typeof(byear_ranges)

## create birthyear_bins
  byear_bins <- paste0("groupingBIN",1:4)
  
## combine [birthyear_bins] and  [birthyear_ranges]
  #birthyear <- cbind(byear_bins, byear_ranges)
  birthyear <- list()
  for (i in 1:length(byear_bins)) {
   birthyear[[byear_bins[i]]] <- byear_ranges[i]  # a list
  } #for: i
  
  print(birthyear)
  
# Set working directory to where your LDSC log files are
setwd(file.path(gc_SI_dir, "/ever_tobacco_user_birthyear"))

# List LDSC log files matching the naming pattern
log_files <- list.files(pattern = "4_rg_raw.log$")
log_files

# Function to extract rg summary from one log file
extract_rg_summary <- function(file) {
  lines <- readLines(file)
  start <- grep("Summary of Genetic Correlation Results", lines)
  if (length(start) == 0) return(NULL)
  data_lines <- lines[(start + 1):length(lines)]
  data_lines <- data_lines[data_lines != ""]
  
  tryCatch({
    df <- read.table(text = data_lines, fill = TRUE, stringsAsFactors = FALSE)
    colnames(df) <- c("p1", "p2", "rg", "se", "z", "p",
                      "h2_obs", "h2_obs_se", "h2_int", "h2_int_se",
                      "gcov_int", "gcov_int_se")
    
    # Extract bin names in format "groupingBIN1", "groupingBIN2", "groupingBIN3", "groupingBIN4"
    extract_bin <- function(path) sub(".*_groupingBIN(\\d+).*", "groupingBIN\\1", path)
    df$Bin1 <- sapply(df$p1, extract_bin)
    df$Bin2 <- sapply(df$p2, extract_bin)
    
    return(df)
  }, error = function(e) NULL)
}


# Extract and combine all data
all_data <- map_df(log_files, extract_rg_summary)
all_data 
           
# Define full bin range
bins <- paste0("groupingBIN", 1:4)

# Filter, pivot, and clean up matrix
# rg
rg_matrix_birthyear <- all_data %>%
                        filter(Bin1 %in% bins, Bin2 %in% bins) %>%
                        distinct(Bin1, Bin2, .keep_all = TRUE) %>%
                        select(Bin1, Bin2, rg) %>%
                        pivot_wider(names_from = Bin2, values_from = rg) %>%
                        column_to_rownames("Bin1") %>%
                        select(any_of(bins)) %>%
                        .[intersect(rownames(.), bins), ]
# se
rg_se_matrix_birthyear <- all_data %>%
                        filter(Bin1 %in% bins, Bin2 %in% bins) %>%
                        distinct(Bin1, Bin2, .keep_all = TRUE) %>%
                        select(Bin1, Bin2, se) %>%
                        pivot_wider(names_from = Bin2, values_from = se) %>%
                        column_to_rownames("Bin1") %>%
                        select(any_of(bins)) %>%
                        .[intersect(rownames(.), bins), ]


# Print the matrix
cat("\n[rg_matrix_birthyear]\n")
print(rg_matrix_birthyear)
cat("\n[rg_se_matrix_birthyear]\n")
print(rg_se_matrix_birthyear)

# Reorder rows
bins_or <- paste0("groupingBIN", 1:4)
rg_matrix_birthyear_ro <- rg_matrix_birthyear[bins_or, ]
rg_se_matrix_birthyear_ro <- rg_se_matrix_birthyear[bins_or, ]

# Print the matrix
cat("\n[rg_matrix_birthyear_ro]\n")
print(rg_matrix_birthyear_ro)
cat("\n[rg_se_matrix_birthyear_ro]\n")
print(rg_se_matrix_birthyear_ro)


# ---------------------------------------------------------
# 1. Extract values
# ---------------------------------------------------------
## heritability
h2_full   <- heritability_full$h2[heritability_full$Sample == "Full"]
h2_female <- heritability_gender$h2[heritability_gender$Sample == "Female"]
h2_male   <- heritability_gender$h2[heritability_gender$Sample == "Male"]
h2_bin1   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN1"]
h2_bin2   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN2"]
h2_bin3   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN3"]
h2_bin4   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN4"]

## genetic correlation
rg_gender_val <- as.numeric(rg_matrix_gender["Female", "Male"])
rg_gender_se  <- as.numeric(rg_se_matrix_gender["Female", "Male"])
rg_bin12_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN1", "groupingBIN2"])
rg_bin12_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN1", "groupingBIN2"])
rg_bin13_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN1", "groupingBIN3"])
rg_bin13_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN1", "groupingBIN3"])
rg_bin14_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN1", "groupingBIN4"])
rg_bin14_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN1", "groupingBIN4"])
rg_bin23_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN2", "groupingBIN3"])
rg_bin23_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN2", "groupingBIN3"])
rg_bin24_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN2", "groupingBIN4"])
rg_bin24_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN2", "groupingBIN4"])
rg_bin34_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN3", "groupingBIN4"])
rg_bin34_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN3", "groupingBIN4"])

# define column order
groups_all <- c("Full", "Female", "Male", "groupingBIN1", "groupingBIN2", "groupingBIN3", "groupingBIN4")
n <- length(groups_all)

# ---------------------------------------------------------
# 2. Build triangle data
# ---------------------------------------------------------
data_tri <- data.frame(
  Row = c("groupingBIN4",
          "groupingBIN3", "groupingBIN3", 
          "groupingBIN2", "groupingBIN2", "groupingBIN2", "groupingBIN2", 
          "groupingBIN1", "groupingBIN1", "groupingBIN1", "groupingBIN1", "groupingBIN1", "groupingBIN1", "groupingBIN1"),
  Col = c(                                                                          "groupingBIN4",
                                                                    "groupingBIN3", "groupingBIN4",
                            "Male",                 "groupingBIN2", "groupingBIN3", "groupingBIN4",
          "Full", "Female", "Male", "groupingBIN1", "groupingBIN2", "groupingBIN3", "groupingBIN4"),
  value = c(
                                                                            h2_bin4,
                                                              h2_bin3,      rg_bin34_val,
                        h2_male,                h2_bin2,      rg_bin23_val, rg_bin24_val,
    h2_full, h2_female, rg_gender_val, h2_bin1, rg_bin12_val, rg_bin13_val, rg_bin14_val)
)

#cat("\n[data_tri]:\n")
#data_tri

# ---------------------------------------------------------
# 3. Add SE
# ---------------------------------------------------------
# add SE columns
data_tri$se <- NA

# assign SE only to rg cells
data_tri$se[data_tri$value == rg_gender_val] <- rg_gender_se
data_tri$se[data_tri$value == rg_bin12_val]  <- rg_bin12_se
data_tri$se[data_tri$value == rg_bin13_val]  <- rg_bin13_se
data_tri$se[data_tri$value == rg_bin14_val]  <- rg_bin14_se
data_tri$se[data_tri$value == rg_bin23_val]  <- rg_bin23_se
data_tri$se[data_tri$value == rg_bin24_val]  <- rg_bin24_se
data_tri$se[data_tri$value == rg_bin34_val]  <- rg_bin34_se

# identify diagonal (heritability) cells
data_tri$is_h2 <- is.na(data_tri$se)  # h² cells have no SE in your construction
cat("\n[data_tri]:\n")
data_tri

# ---------------------------------------------------------
# 4. Significance test for rg = 1
# ---------------------------------------------------------
# compute stars for all rg values
data_tri$star <- ""
rg_rows <- !is.na(data_tri$se)

data_tri <- data_tri %>%
                mutate(
                    n_tests = ifelse(!is_h2&(Col=="Male"), 1, sum(rg_rows)-1),
                    z_vals = ifelse(!is_h2, (value - 1) / se, NA_real_),
                    p_vals = ifelse(!is_h2, 2 * pnorm(-abs(z_vals)), NA_real_),
                    star = ifelse(!is_h2&(Col!="Male"), 
                                  ifelse(p_vals < 0.05/n_tests, "**",
                                        ifelse(p_vals < 0.05, "*")
                                        ) #ifelse
                                  , ifelse(!is_h2&(Col=="Male"), "*"
                                           , NA_real_)),
                    label = ifelse(is.na(value), "",
                                   ifelse(is_h2,
                                          sprintf("h²=%.3f", value),
                                          sprintf("%.3f%s", value, star)
                                         )  #ifelse
                                  ) #ifelse
                ) #mutate

cat("\n[data_tri]:\n")
data_tri


# ---------------------------------------------------------
# 5. Plot
# ---------------------------------------------------------
p_h2_rg_stair <- ggplot(data_tri, aes(Col, Row)) +
  geom_tile(aes(fill = ifelse(!is.na(se), value, NA)), color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.2) +
  geom_vline(xintercept = c(1.5, 3.5), color = "white", linewidth = 1.9) +
  scale_fill_gradient2(
    low = "#A020F0",   # bright purple
    mid = "white",
    high = "#FF0000",  # bright red
    midpoint = 0.5,
    limits = c(0, 1),
    na.value = "white"   # important
  ) +
  scale_x_discrete(limits = groups_all, drop = T) +
  scale_y_discrete(drop = T, position = "right",
                   labels = function(x) {
                     ifelse(x == "groupingBIN1", "Silent Generation",
                     ifelse(x == "groupingBIN2", "Baby Boomers",
                     ifelse(x == "groupingBIN3", "Generation X", "")))
                   }) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
#    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
                axis.text.y.right = element_text(hjust = -0.15),
      legend.position = "none"
  )

p_h2_rg_stair

ggsave(paste0(working_dir, "/output/Figures/h2_rg_stair.png"), plot=p_h2_rg_stair, width = 10, height = 6, dpi = 500)

# ---------------------------------------------------------
# 1. Extract values
# ---------------------------------------------------------
## heritability
h2_full   <- heritability_full$h2[heritability_full$Sample == "Full"]
h2_female <- heritability_gender$h2[heritability_gender$Sample == "Female"]
h2_male   <- heritability_gender$h2[heritability_gender$Sample == "Male"]
h2_bin1   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN1"]
h2_bin2   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN2"]
h2_bin3   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN3"]
h2_bin4   <- heritability_generation$h2[heritability_generation$Sample == "groupingBIN4"]

## genetic correlation
rg_gender_val <- as.numeric(rg_matrix_gender["Female", "Male"])
rg_gender_se  <- as.numeric(rg_se_matrix_gender["Female", "Male"])
rg_bin12_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN1", "groupingBIN2"])
rg_bin12_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN1", "groupingBIN2"])
rg_bin13_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN1", "groupingBIN3"])
rg_bin13_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN1", "groupingBIN3"])
rg_bin14_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN1", "groupingBIN4"])
rg_bin14_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN1", "groupingBIN4"])
rg_bin23_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN2", "groupingBIN3"])
rg_bin23_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN2", "groupingBIN3"])
rg_bin24_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN2", "groupingBIN4"])
rg_bin24_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN2", "groupingBIN4"])
rg_bin34_val <- as.numeric(rg_matrix_birthyear_ro["groupingBIN3", "groupingBIN4"])
rg_bin34_se  <- as.numeric(rg_se_matrix_birthyear_ro["groupingBIN3", "groupingBIN4"])

# ---------------------------------------------------------
# 2. Build triangle data
# ---------------------------------------------------------
data_tri <- data.frame(
  Row = c("groupingBIN4",
          "groupingBIN3","groupingBIN3",
          "groupingBIN2","groupingBIN2","groupingBIN2","groupingBIN2",
          "groupingBIN1","groupingBIN1","groupingBIN1",
          "groupingBIN1","groupingBIN1","groupingBIN1","groupingBIN1"),
  Col = c("groupingBIN4",
          "groupingBIN3","groupingBIN4",
          "Male","groupingBIN2","groupingBIN3","groupingBIN4",
          "Full","Female","Male",
          "groupingBIN1","groupingBIN2","groupingBIN3","groupingBIN4"),
  value = c(h2_bin4,
            h2_bin3,rg_bin34_val,
            h2_male,h2_bin2,rg_bin23_val,rg_bin24_val,
            h2_full,h2_female,rg_gender_val,
            h2_bin1,rg_bin12_val,rg_bin13_val,rg_bin14_val)
)

# ---------------------------------------------------------
# 3. Recode BIN → Gen (match bottom heatmap)
# ---------------------------------------------------------

data_tri <- data_tri %>%
  mutate(
    Row = recode(Row,
                 groupingBIN1="Gen1",
                 groupingBIN2="Gen2",
                 groupingBIN3="Gen3",
                 groupingBIN4="Gen4"),
    Col = recode(Col,
                 groupingBIN1="Gen1",
                 groupingBIN2="Gen2",
                 groupingBIN3="Gen3",
                 groupingBIN4="Gen4")
  )

# ---------------------------------------------------------
# 4. Add SE
# ---------------------------------------------------------

data_tri$se <- NA
data_tri$se[data_tri$value == rg_gender_val] <- rg_gender_se
data_tri$se[data_tri$value == rg_bin12_val]  <- rg_bin12_se
data_tri$se[data_tri$value == rg_bin13_val]  <- rg_bin13_se
data_tri$se[data_tri$value == rg_bin14_val]  <- rg_bin14_se
data_tri$se[data_tri$value == rg_bin23_val]  <- rg_bin23_se
data_tri$se[data_tri$value == rg_bin24_val]  <- rg_bin24_se
data_tri$se[data_tri$value == rg_bin34_val]  <- rg_bin34_se

data_tri$is_h2 <- is.na(data_tri$se)

# ---------------------------------------------------------
# 5. Significance test for rg = 1
# ---------------------------------------------------------
# compute stars for all rg values
data_tri$star <- ""
rg_rows <- !is.na(data_tri$se)

data_tri <- data_tri %>%
                mutate(
                    n_tests = ifelse(!is_h2&(Col=="Male"), 1, sum(rg_rows)-1),
                    z_vals = ifelse(!is_h2, (value - 1) / se, NA_real_),
                    p_vals = ifelse(!is_h2, 2 * pnorm(-abs(z_vals)), NA_real_),
                    star = ifelse(!is_h2&(Col!="Male"), 
                                  ifelse(p_vals < 0.05/n_tests, "**",
                                        ifelse(p_vals < 0.05, "*")
                                        ) #ifelse
                                  , ifelse(!is_h2&(Col=="Male"), "*"
                                           , NA_real_)),
                    label = ifelse(is.na(value), "",
                                   ifelse(is_h2,
                                          sprintf("h²=%.3f", value),
                                          sprintf("%.3f%s", value, star)
                                         )  #ifelse
                                  ) #ifelse
                ) #mutate

cat("\n[data_tri]:\n")
data_tri

# ---------------------------------------------------------
# 6. Align x-axis with bottom heatmap + spacer
# ---------------------------------------------------------

shared_x <- c(
  "Full", "Spacer", "Female", "Male", "F vs M",
  "Gen1", "Gen2", "Gen3", "Gen4", "Gen Joint"
)

data_tri$Col <- factor(data_tri$Col, levels = shared_x)
data_tri$Row <- factor(data_tri$Row, levels = c("Gen1","Gen2","Gen3","Gen4"))

# Add gray spacer column
data_spacer_tri <- data.frame(
  Row = factor(c("Gen1", "Gen2", "Gen3"), 
               levels = c("Gen1","Gen2","Gen3","Gen4")),
  Col = factor("Spacer", levels = shared_x),
  value = NA_real_,
  se = NA_real_,
  is_h2 = TRUE,
  star = "",
  label = ""
)

data_tri2 <- bind_rows(data_tri, data_spacer_tri)


# Female label only
label_female <- data.frame(
  Row = factor("Gen1",
               levels = c("Gen1","Gen2","Gen3","Gen4")),
  Col = factor("F vs M",
               levels = shared_x),
  label = "Female"
)
# ---------------------------------------------------------
# 7. Plot
# ---------------------------------------------------------
p_h2_rg_stair_new <- ggplot(data_tri2, aes(Col, Row)) +
  geom_tile(aes(fill = ifelse(!is.na(se), value, NA)), color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.2) +
  geom_text(data = label_female, aes(Col, Row, label = label), inherit.aes = FALSE,
            size = 2.85, color = "grey30") +
  scale_fill_gradient2(
    low = "#A020F0",   # bright purple
    mid = "white",
    high = "#FF0000",  # bright red
    midpoint = 0.5,
    limits = c(0, 1),
    na.value = "white"   # important
  ) +
  scale_x_discrete(limits = shared_x, drop = FALSE,
                  labels = function(x) {
                                          x[x == "Spacer"] <- ""
                                          x
                                        }) +
  scale_y_discrete(drop = TRUE, position = "right",
                   labels = function(x) {
                     ifelse(x == "Gen1", "Silent Generation",
                     ifelse(x == "Gen2", "Baby Boomers",
                     ifelse(x == "Gen3", "Generation X", "")))
                   }) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y.right = element_text(hjust = 0, color="grey30", size=7.55, margin = margin(l = -65)),
    axis.text.y.left = element_blank(),
      legend.position = "none"
  )

p_h2_rg_stair_new

ggsave(paste0(working_dir, "/output/Figures/h2_rg_stair_new.png"), plot=p_h2_rg_stair_new, width = 10, height = 6, dpi = 500)

final_plot <- p_h2_rg_stair / gc_traits +
                plot_layout(heights = c(0.85, 4.5)) &  #c(1.2, 3.8))  # adjust ratio here
                theme(plot.margin = margin(0, 0, 7, 0))

final_plot

ggsave(paste0(working_dir, "/output/Figures/h2_rg_stair_gc_traits.png"), 
       plot=final_plot, width = 10, height = 9, dpi = 500)

final_plot <- p_h2_rg_stair_new / gc_traits_new +
                plot_layout(heights = c(0.75, 4.5)) &  #c(1.2, 3.8))  # adjust ratio here
                theme(plot.margin = margin(0, 0, 7, 0))

final_plot

ggsave(paste0(working_dir, "/output/Figures/h2_rg_stair_gc_traits_new.png"), 
       plot=final_plot, width = 10, height = 9, dpi = 500)


