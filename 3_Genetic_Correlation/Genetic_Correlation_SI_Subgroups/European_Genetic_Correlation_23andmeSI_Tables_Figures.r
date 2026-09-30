R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()  

#Define directories for "EUROPEAN"
# 1) Define Ancestry
ANCESTRY="european"

# 1) Define directories
cat("Set directory\n")
root_dir=paste0("/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC")
results_dir=paste0(root_dir,"/Genetic_Correlation/raw/SI")
results_dir_tables=paste0(results_dir,"/Tables")
results_dir_figures=paste0(results_dir,"/Figures")
cat("\n-------------------------------------------------------------------------------\n")

#Current directory
cat("Current directory\n")
setwd(root_dir)
cat(getwd())
cat("\n-------------------------------------------------------------------------------\n")

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
cat("\n*Load [knitr]\n") 
library(knitr)
cat("\n")
cat("\n*Load [tidyverse]\n") 
library(tidyverse)
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
cat("\n-----------------------------------------------------------------------------------\n")


#-------------------------------- Gender --------------------------------------#
# Set working directory to where your LDSC log files are
setwd(file.path(results_dir, "/ever_tobacco_user_gender"))

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

#Test for perfect genetic correlation (H0: rg=1)
rg_hat <- as.numeric(rg_matrix_gender)
rg_se <- as.numeric(rg_se_matrix_gender)
print(c(rg_hat, rg_se))

test_rg <- function(eq, rg, se) {
    z <- (rg - eq)/se
    p <- 2*pnorm(-abs(z))  # pnorm() is a built-in function used to calculate the CDF of the normal dist.
    return(list(z = z, p_value = p))
} #test_rg
cat("\n")
cat("[H0: rg=0]")
test_rg(0, rg_hat, rg_se)
cat("[H0: rg=1]")
test_rg(1, rg_hat, rg_se)

#Get sig stars
get_sig_star <- function(p, num_test, bonf = FALSE) {
    if (bonf) {
        if (p < 0.05/num_test) return("**")
        if (p < 0.05) return("*")
    } else {
        if (p < 0.05) return("*")
    } # if-else
    return("")
} # get_sig_star

#Test H0: rg=1
test <- test_rg(1, rg_hat, rg_se)
#Get asterisk (no Bonferroni for sex)
star <- get_sig_star(test$p_value, 1, bonf=F)
cat("star")
star

#------------------- Create [Genetic Correlation] Tables ----------------------#
## Gender sample

# Convert matrix to data frame
rg_matrix_gender_clean <- as.data.frame(rg_matrix_gender)

# Replace NA with empty string
rg_matrix_gender_clean[is.na(rg_matrix_gender_clean)] <- ""

# View results
print(rg_matrix_gender_clean)

# Create LaTeX table object
xtab <- xtable(
  rg_matrix_gender_clean,
  caption = "Genetic Correlation (Gender)",
  label = "tab:rg_gender"
)

print(
  xtab,
  file = paste0(results_dir_tables, "/rg_gender.tex"),
  include.rownames = TRUE,
  caption.placement = "top",
  sanitize.text.function = identity,
  floating = TRUE,
  booktabs = TRUE
)
sink()

#------------------- Create [Genetic Correlation] Plots------------------------#

# Melt into long format
df_long <- melt(as.matrix(rg_matrix_gender), na.rm = TRUE)

# Ensure value is numeric
df_long$value <- as.numeric(df_long$value)

# Plot
ggplot(df_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(round(value, 4),star)), size = 6) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0, 1), name = "rg") +
  theme_minimal() +
  labs(title = "Genetic Correlation",
       subtitle = "By Sex",
       x = "",
       y = "") +
  theme_minimal(base_size = 18) +  # Minimalistic theme
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/rg_gender.png"), width = 8, height = 6, dpi = 500)

#------------------------------------ Region ---------------------------------------#
# Set working directory to where your LDSC log files are
setwd(file.path(results_dir, "/ever_tobacco_user_region"))

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

#Test for perfect genetic correlation (H0: rg=1)

# Convert to numeric matrices
rg_hat_mat <- as.matrix(rg_matrix_region)
rg_se_mat <- as.matrix(rg_se_matrix_region)
 cat("[rg_hat_mat]\n")
 print(rg_hat_mat)
 cat("[rg_se_mat]\n")
 print(rg_se_mat)
 cat("\n")

# Identify valid (non-NA) entries
valid_idx <- which(!is.na(rg_hat_mat), arr.ind = TRUE)
valid_idx

# Extract vectors
rg_hat <- as.numeric(rg_hat_mat[valid_idx])
rg_se <- as.numeric(rg_se_mat[valid_idx])
 cat("[rg_hat]\n")
 print(rg_hat)
 cat("[rg_se]\n")
 print(rg_se)
 cat("\n")

cat("length: [rg_hat]")
length(rg_hat)

# Test function
test_rg <- function(eq, rg, se) {
    z <- (rg - eq)/se
    p <- 2*pnorm(-abs(z))  # pnorm() is a built-in function used to calculate the CDF of the normal dist.
    return(list(z = z, p_value = p))
} #test_rg
cat("\n")
cat("[H0: rg=0]")
test_rg(0, rg_hat, rg_se)
cat("[H0: rg=1]")
test_rg(1, rg_hat, rg_se)

# Test H0: rg=1
test <- test_rg(1, rg_hat, rg_se)
p_vals <- test$p_value
p_vals

# Number of tests
num_tests <- length(p_vals)
num_tests

#Get sig stars
get_sig_star <- function(p, num_test, bonf = FALSE) {
    if (bonf) {
        if (p < 0.05/num_test) return("**")
        if (p < 0.05) return("*")
    } else {
        if (p < 0.05) return("*")
    } # if-else
    return("")
} # get_sig_star

#Get asterisk (no Bonferroni for sex)
star <- sapply(test$p_value, get_sig_star, num_test = num_tests, bonf = TRUE) 
cat("[star]")
star

#Put stars back into a matrix
# Initializa star matrix
star_mat <- matrix(NA, nrow = nrow(rg_hat_mat), ncol = ncol(rg_hat_mat))
rownames(star_mat) <- rownames(rg_hat_mat)
colnames(star_mat) <- colnames(rg_hat_mat)
star_mat
# Fill stars
for (i in seq_along(star)) {
  star_mat[valid_idx[i, 1], valid_idx[i, 2]] <- star[i]
} #for: i
cat("[star_mat]\n")
star_mat

# Create label matrix (rg+star): This is what I will feed to [geom_text()]
label_mat <- matrix("", nrow = nrow(rg_hat_mat), ncol = ncol(rg_hat_mat))
rownames(label_mat) <- rownames(rg_hat_mat)
colnames(label_mat) <- colnames(rg_hat_mat)
label_mat
# Fill in label matrix
label_mat[!is.na(rg_hat_mat)] <- paste0(round(as.numeric(rg_hat_mat[!is.na(rg_hat_mat)]), 4), 
                                        star_mat[!is.na(rg_hat_mat)])
cat("[label_mat]\n")
label_mat

#------------------- Create [Genetic Correlation] Tables ----------------------#
## REGION sample

# reverse the row order to align with the graph.
rg_matrix_region_rr <- rg_matrix_region[nrow(rg_matrix_region):1, ]

# Print the row reversed matrix
print(rg_matrix_region_rr)

# Convert matrix to data frame and remove NAs by replacing with ""
rg_matrix_region_rr_clean <- as.data.frame(rg_matrix_region_rr)
rg_matrix_region_rr_clean[is.na(rg_matrix_region_rr_clean)] <- ""

# Create LaTeX table object
xtab_region <- xtable(
  rg_matrix_region_rr_clean,
  caption = "Genetic Correlation (Region)",
  label = "tab:rg_region"
)

print(
  xtab_region,
  file = paste0(results_dir_tables, "/rg_region.tex"),
  include.rownames = TRUE,
  caption.placement = "top",
  sanitize.text.function = identity,
  floating = TRUE,
  booktabs = TRUE
)

#------------------- Create [Genetic Correlation] Plots------------------------#

# Melt [rg] matrix (only non-NA pairs)
df_long_temp <- melt(as.matrix(rg_matrix_region), na.rm = TRUE)
 #df_long_temp
# Melt [label] matrix WITHOUT dropping anything, then keep only matching pairs
df_label <- melt(label_mat, na.rm = TRUE)
colnames(df_label) <- c("Var1", "Var2", "label")
 #df_label
# Join [df_long] and [df_label]
df_long <- df_long_temp %>%
                left_join(df_label, by = c("Var1", "Var2"))

 #df_long

# Ensure value is numeric
df_long$value <- as.numeric(df_long$value)
df_long

# Plot
ggplot(df_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 6) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0, 1), name = "rg") +
  theme_minimal() +
  labs(title = "Genetic Correlation",
       subtitle = "By Region",
       x = "",
       y = "") +
  theme_minimal(base_size = 18) +  # Minimalistic theme
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/rg_region.png"), width = 8, height = 6, dpi = 300)


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
setwd(file.path(results_dir, "/ever_tobacco_user_birthyear"))

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


#Test for perfect genetic correlation (H0: rg=1)

# Convert to numeric matrices
rg_hat_mat <- as.matrix(rg_matrix_birthyear_ro)
rg_se_mat <- as.matrix(rg_se_matrix_birthyear_ro)
 cat("[rg_hat_mat]\n")
 print(rg_hat_mat)
 cat("[rg_se_mat]\n")
 print(rg_se_mat)
 cat("\n")

# Identify valid (non-NA) entries
valid_idx <- which(!is.na(rg_hat_mat), arr.ind = TRUE)
valid_idx

# Extract vectors
rg_hat <- as.numeric(rg_hat_mat[valid_idx])
rg_se <- as.numeric(rg_se_mat[valid_idx])
 cat("[rg_hat]\n")
 print(rg_hat)
 cat("[rg_se]\n")
 print(rg_se)
 cat("\n")

cat("length: [rg_hat]")
length(rg_hat)

# Test function
test_rg <- function(eq, rg, se) {
    z <- (rg - eq)/se
    p <- 2*pnorm(-abs(z))  # pnorm() is a built-in function used to calculate the CDF of the normal dist.
    return(list(z = z, p_value = p))
} #test_rg
cat("\n")
cat("[H0: rg=0]")
test_rg(0, rg_hat, rg_se)
cat("[H0: rg=1]")
test_rg(1, rg_hat, rg_se)

# Test H0: rg=1
test <- test_rg(1, rg_hat, rg_se)
p_vals <- test$p_value
p_vals

# Number of tests
num_tests <- length(p_vals)
num_tests

#Get sig stars
get_sig_star <- function(p, num_test, bonf = FALSE) {
    if (bonf) {
        if (p < 0.05/num_test) return("**")
        if (p < 0.05) return("*")
    } else {
        if (p < 0.05) return("*")
    } # if-else
    return("")
} # get_sig_star

#Get asterisk (no Bonferroni for sex)
star <- sapply(test$p_value, get_sig_star, num_test = num_tests, bonf = TRUE) 
cat("[star]")
star

#Put stars back into a matrix
# Initializa star matrix
star_mat <- matrix(NA, nrow = nrow(rg_hat_mat), ncol = ncol(rg_hat_mat))
rownames(star_mat) <- rownames(rg_hat_mat)
colnames(star_mat) <- colnames(rg_hat_mat)
star_mat
# Fill stars
for (i in seq_along(star)) {
  star_mat[valid_idx[i, 1], valid_idx[i, 2]] <- star[i]
} #for: i
cat("[star_mat]\n")
star_mat

# Create label matrix (rg+star): This is what I will feed to [geom_text()]
label_mat <- matrix("", nrow = nrow(rg_hat_mat), ncol = ncol(rg_hat_mat))
rownames(label_mat) <- rownames(rg_hat_mat)
colnames(label_mat) <- colnames(rg_hat_mat)
label_mat
# Fill in label matrix
label_mat[!is.na(rg_hat_mat)] <- paste0(round(as.numeric(rg_hat_mat[!is.na(rg_hat_mat)]), 3), 
                                        star_mat[!is.na(rg_hat_mat)])
cat("[label_mat]\n")
label_mat


#------------------- Create [Genetic Correlation] Tables ----------------------#
## Birthyear sample

# Reorder rows
bins_or2 <- paste0("bin", 18:4)
rg_matrix_birthyear_ro2 <- rg_matrix_birthyear[bins_or2, ]

# Print the matrix
print(rg_matrix_birthyear_ro2)

# Convert matrix to data frame and remove NAs by replacing with ""
rg_matrix_birthyear_ro2_clean <- as.data.frame(rg_matrix_birthyear_ro2)
rg_matrix_birthyear_ro2_clean[is.na(rg_matrix_birthyear_ro2_clean)] <- ""

# Create LaTeX table object
xtab_birthyear <- xtable(
  rg_matrix_birthyear_ro2_clean,
  caption = "Genetic Correlation (Birthyear Bin)",
  label = "tab:rg_birthyear"
)

print(
  xtab_birthyear,
  file = paste0(results_dir_tables, "/rg_birthyear.tex"),
  include.rownames = TRUE,
  caption.placement = "top",
  sanitize.text.function = identity,
  floating = TRUE,
  booktabs = TRUE
)


#------------------- Create [Genetic Correlation] Plots------------------------#

# Melt [rg] matrix (only non-NA pairs)
df_long_temp <- melt(as.matrix(rg_matrix_birthyear_ro), na.rm = TRUE)
# Melt [label] matrix WITHOUT dropping anything, then keep only matching pairs
df_label <- melt(label_mat, na.rm = TRUE)
colnames(df_label) <- c("Var1", "Var2", "label")
# Join [df_long] and [df_label]
df_long <- df_long_temp %>%
                left_join(df_label, by = c("Var1", "Var2"))

# Ensure value is numeric
df_long$value <- as.numeric(df_long$value)
df_long 

# [bin] to [birthyear interval]
 # Create a named vector for bin-to-birthyear range mapping
 bin_labels <- unlist(birthyear)  # Convert list to named character vector
 names(bin_labels) <- names(birthyear)  # Just in case  # Now bin_labels["bin4"] returns "1924-1928"

 #Map axis labels in df_long
 df_long$Var1_label <- bin_labels[as.character(df_long$Var1)]
 df_long$Var2_label <- bin_labels[as.character(df_long$Var2)]


# Plot
ggplot(df_long, aes(x = Var2_label, y = Var1_label, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0, 1), name = "rg") +
  theme_minimal() +
  labs(title = "Genetic Correlation",
       subtitle = "By Birth Year",
       x = "",
       y = "") +
  theme_minimal(base_size = 18) +  # Minimalistic theme
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )


 # Save the plot
 ggsave(paste0(results_dir_figures,"/rg_birthyear.png"), width = 8, height = 6, dpi = 500)

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
setwd(file.path(results_dir, "/ever_tobacco_user_birthyear"))

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


#Test for perfect genetic correlation (H0: rg=1)

# Convert to numeric matrices
rg_hat_mat <- as.matrix(rg_matrix_birthyear_ro)
rg_se_mat <- as.matrix(rg_se_matrix_birthyear_ro)
 cat("[rg_hat_mat]\n")
 print(rg_hat_mat)
 cat("[rg_se_mat]\n")
 print(rg_se_mat)
 cat("\n")

# Identify valid (non-NA) entries
valid_idx <- which(!is.na(rg_hat_mat), arr.ind = TRUE)
valid_idx

# Extract vectors
rg_hat <- as.numeric(rg_hat_mat[valid_idx])
rg_se <- as.numeric(rg_se_mat[valid_idx])
 cat("[rg_hat]\n")
 print(rg_hat)
 cat("[rg_se]\n")
 print(rg_se)
 cat("\n")

cat("length: [rg_hat]")
length(rg_hat)

# Test function
test_rg <- function(eq, rg, se) {
    z <- (rg - eq)/se
    p <- 2*pnorm(-abs(z))  # pnorm() is a built-in function used to calculate the CDF of the normal dist.
    return(list(z = z, p_value = p))
} #test_rg
cat("\n")
cat("[H0: rg=0]")
test_rg(0, rg_hat, rg_se)
cat("[H0: rg=1]")
test_rg(1, rg_hat, rg_se)

# Test H0: rg=1
test <- test_rg(1, rg_hat, rg_se)
p_vals <- test$p_value
p_vals

# Number of tests
num_tests <- length(p_vals)
num_tests

#Get sig stars
get_sig_star <- function(p, num_test, bonf = FALSE) {
    if (bonf) {
        if (p < 0.05/num_test) return("**")
        if (p < 0.05) return("*")
    } else {
        if (p < 0.05) return("*")
    } # if-else
    return("")
} # get_sig_star

#Get asterisk (no Bonferroni for sex)
star <- sapply(test$p_value, get_sig_star, num_test = num_tests, bonf = TRUE) 
cat("[star]")
star

#Put stars back into a matrix
# Initializa star matrix
star_mat <- matrix(NA, nrow = nrow(rg_hat_mat), ncol = ncol(rg_hat_mat))
rownames(star_mat) <- rownames(rg_hat_mat)
colnames(star_mat) <- colnames(rg_hat_mat)
star_mat
# Fill stars
for (i in seq_along(star)) {
  star_mat[valid_idx[i, 1], valid_idx[i, 2]] <- star[i]
} #for: i
cat("[star_mat]\n")
star_mat

# Create label matrix (rg+star): This is what I will feed to [geom_text()]
label_mat <- matrix("", nrow = nrow(rg_hat_mat), ncol = ncol(rg_hat_mat))
rownames(label_mat) <- rownames(rg_hat_mat)
colnames(label_mat) <- colnames(rg_hat_mat)
label_mat
# Fill in label matrix
label_mat[!is.na(rg_hat_mat)] <- paste0(round(as.numeric(rg_hat_mat[!is.na(rg_hat_mat)]), 4), 
                                        star_mat[!is.na(rg_hat_mat)])
cat("[label_mat]\n")
label_mat


#------------------- Create [Genetic Correlation] Tables ----------------------#
## Birthyear sample

# Reorder rows
bins_or2 <- paste0("groupingBIN", 3:1)
rg_matrix_birthyear_ro2 <- rg_matrix_birthyear[bins_or2, ]

# Print the matrix
print(rg_matrix_birthyear_ro2)

# Convert matrix to data frame and remove NAs by replacing with ""
rg_matrix_birthyear_ro2_clean <- as.data.frame(rg_matrix_birthyear_ro2)
rg_matrix_birthyear_ro2_clean[is.na(rg_matrix_birthyear_ro2_clean)] <- ""

# Create LaTeX table object
xtab_birthyear <- xtable(
  rg_matrix_birthyear_ro2_clean,
  caption = "Genetic Correlation (Birthyear Bin)",
  label = "tab:rg_birthyear"
)

print(
  xtab_birthyear,
  file = paste0(results_dir_tables, "/rg_birthyear_groupingBINs.tex"),
  include.rownames = TRUE,
  caption.placement = "top",
  sanitize.text.function = identity,
  floating = TRUE,
  booktabs = TRUE
)


#------------------- Create [Genetic Correlation] Plots------------------------#

# Melt [rg] matrix (only non-NA pairs)
df_long_temp <- melt(as.matrix(rg_matrix_birthyear_ro), na.rm = TRUE)
# Melt [label] matrix WITHOUT dropping anything, then keep only matching pairs
df_label <- melt(label_mat, na.rm = TRUE)
colnames(df_label) <- c("Var1", "Var2", "label")
# Join [df_long] and [df_label]
df_long <- df_long_temp %>%
                left_join(df_label, by = c("Var1", "Var2"))

# Ensure value is numeric
df_long$value <- as.numeric(df_long$value)
df_long

# [bin] to [birthyear interval]
 # Create a named vector for bin-to-birthyear range mapping
 bin_labels <- unlist(birthyear)  # Convert list to named character vector
 names(bin_labels) <- names(birthyear)  # Just in case  # Now bin_labels["bin4"] returns "1924-1928"

 #Map axis labels in df_long
 df_long$Var1_label <- bin_labels[as.character(df_long$Var1)]
 df_long$Var2_label <- bin_labels[as.character(df_long$Var2)]


# Plot
ggplot(df_long, aes(x = Var2_label, y = Var1_label, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 6) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0, 1), name = "rg") +
  theme_minimal() +
  labs(title = "Genetic Correlation",
       subtitle = "By Birth Year",
       x = "",
       y = "") +
  theme_minimal(base_size = 18) +  # Minimalistic theme
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )


 # Save the plot
 ggsave(paste0(results_dir_figures,"/rg_birthyear_groupingBINs.png"), width = 8, height = 6, dpi = 500)




