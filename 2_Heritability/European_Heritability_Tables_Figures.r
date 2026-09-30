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
results_dir=paste0(root_dir,"/Heritability/raw")
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
#Libraries useful for creating latex tables (.tex file)
cat("\n*Load [stargazer]\n") #For regression tables
# install.packages("stargazer")
library(stargazer)
cat("\n")
cat("\n*Load [kableExtra]\n") #For flexible formatting
# install.packages("kableExtra")
#library(kableExtra)
cat("\n")
cat("\n*Load [xtable]\n") #For a simple approach
library(xtable)
cat("\n")
cat("\n*Load [knitr]\n") #For super column
#library(knitr)
cat("\n")
cat("\n*Load [purrr]\n") 
library(purrr)
cat("\n")
cat("\n*Load [ggplot2]\n") 
library(ggplot2)
cat("\n")
cat("\n-----------------------------------------------------------------------------------\n")


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
log_files <- file.path(results_dir, paste0(subdirs, "_h2_raw.log"))


# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# Rename Sample
heritability_df$Sample[heritability_df$Sample == "full"] <- "Full"

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# View results
print(heritability_df)

#---------------------- Create [Heritability] Tables --------------------------#
## Full and Gender sample

# Transpose
heritability_full <- as.data.frame(t(heritability_df))

# View results
print(heritability_full)

# Create LaTeX table object
latex_table <- xtable(heritability_full,
                      caption = "Heritability (Full Sample)")

# Optional: customize alignment if needed
# For example, 3 columns with left alignment:
# align = c("l", "l", "l", "l")

# Save LaTeX table to .tex file
sink(file = file.path(results_dir_tables, "h2_full.tex"))
print(latex_table,
      include.rownames = T,
      caption.placement = "top",
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(heritability_full)))  # horizontal lines
sink()


#----------------------- Create [Heritability] Plots---------------------------#

#Create confidence intervals 
heritability_df <- heritability_df %>%
                      mutate(
                             ci_lower = h2 - 1.96 * h2_se,
                             ci_upper = h2 + 1.96 * h2_se
                            )

# set the order of X-axis
heritability_df$Sample <- factor(heritability_df$Sample, levels = c("Full"))

# Draw a plot
ggplot(heritability_df, aes(x = Sample, y = h2)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
#  geom_text(aes(label = round(h2, 4)), vjust = -1.8, hjust = 0.5, size = 6) +
  geom_text(aes(y = ci_upper + 0.005, label = round(h2, 3)), size = 6) +
  labs(
    title = "Observed-Scale Heritability with 95% CI",
    subtitle = "Full Sample",
    x = "",
    y = "Heritability (h²)"
  ) +
  theme_minimal(base_size = 18) +
  ylim(0, max(heritability_df$ci_upper, na.rm = TRUE) + 0.01) +
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/h2_full.png"), width = 8, height = 6, dpi = 500)

# ---------------------------- FULL and GENDER sample -------------------------#

# Target subfolders
subdirs <- c("ever_tobacco_user_female", "ever_tobacco_user_male")

# Construct full log file paths
log_files <- file.path(results_dir, paste0(subdirs, "_h2_raw.log"))


# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# Rename Sample
heritability_df$Sample[heritability_df$Sample == "female"] <- "Female"
heritability_df$Sample[heritability_df$Sample == "male"] <- "Male"

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# View results
print(heritability_df)

#---------------------- Create [Heritability] Tables --------------------------#
## Gender sample

# Transpose
heritability_gender <- as.data.frame(t(heritability_df))

# View results
print(heritability_gender)

# Create LaTeX table object
latex_table <- xtable(heritability_gender,
                      caption = "Heritability (Sex)")

# Optional: customize alignment if needed
# For example, 3 columns with left alignment:
# align = c("l", "l", "l", "l")

# Save LaTeX table to .tex file
sink(file = file.path(results_dir_tables, "h2_gender.tex"))
print(latex_table,
      include.rownames = T,
      caption.placement = "top",
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(heritability_gender)))  # horizontal lines
sink()


# test equality of h2
wald_test <- function(h2_1, se_1, h2_2, se_2) {
  z <- (h2_1 - h2_2) / sqrt(se_1^2 + se_2^2)
  p <- 2 * pnorm(-abs(z))
  return(c(z = z, p_value = p))
}

data <- copy(heritability_gender)

# female vs male
h2_1<-as.numeric(data[2,1])
 h2_1
se_1<-as.numeric(data[3,1])
 se_1
h2_2<-as.numeric(data[2,2])
 h2_2
se_2<-as.numeric(data[3,2])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)

#----------------------- Create [Heritability] Plots---------------------------#

#Create confidence intervals 
heritability_df <- heritability_df %>%
                      mutate(
                             ci_lower = h2 - 1.96 * h2_se,
                             ci_upper = h2 + 1.96 * h2_se
                            )

# set the order of X-axis
heritability_df$Sample <- factor(heritability_df$Sample, levels = c("Female", "Male"))

# Draw a plot
ggplot(heritability_df, aes(x = Sample, y = h2)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
#  geom_text(aes(label = round(h2, 4)), vjust = -1.8, hjust = 0.5, size = 6) +
  geom_text(aes(y = ci_upper + 0.005, label = round(h2, 3)), size = 6) +
  labs(
    title = "Observed-Scale Heritability with 95% CI",
    subtitle = "By Sex",
    x = "",
    y = "Heritability (h²)"
  ) +
  theme_minimal(base_size = 18) +
  ylim(0, max(heritability_df$ci_upper, na.rm = TRUE) + 0.01) +
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/h2_gender.png"), width = 8, height = 6, dpi = 500)

# ------------------------------- REGION sample -------------------------------#

# Target subfolders
subdirs <- paste0("ever_tobacco_user_", c("MW","NE","SE","SW","W"))

# Construct full log file paths
log_files <- file.path(results_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# View results
print(heritability_df)

#---------------------- Create [Heritability] Tables --------------------------#
## Full and Gender sample

# Transpose
heritability_region <- as.data.frame(t(heritability_df))

# View results
print(heritability_region)

# Create LaTeX table object
latex_table <- xtable(heritability_region,
                      caption = "Heritability (Region Sample)")

# Save to a .tex file
sink(file = file.path(results_dir_tables, "h2_region.tex"))
print(latex_table,
      include.rownames = T,
      caption.placement = "top",
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(heritability_region)))
sink()

n_com <- 10
# test equality of h2
wald_test <- function(h2_1, se_1, h2_2, se_2) {
  z <- (h2_1 - h2_2) / sqrt(se_1^2 + se_2^2)
  p <- 2 * pnorm(-abs(z))
  p_bon <- p*n_com
  return(c(z = z, p_value = p, p_value_bon = p_bon))
}

data <- copy(heritability_region)

# MW vs NE
cat("[MW vs NE]")
h2_1<-as.numeric(data[2,1]); h2_1
se_1<-as.numeric(data[3,1]); se_1
h2_2<-as.numeric(data[2,2]); h2_2
se_2<-as.numeric(data[3,2]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# MW vs SE
cat("[MW vs SE]")
h2_1<-as.numeric(data[2,1]); h2_1
se_1<-as.numeric(data[3,1]); se_1
h2_2<-as.numeric(data[2,3]); h2_2
se_2<-as.numeric(data[3,3]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# MW vs SW
cat("[MW vs SW]")
h2_1<-as.numeric(data[2,1]); h2_1
se_1<-as.numeric(data[3,1]); se_1
h2_2<-as.numeric(data[2,4]); h2_2
se_2<-as.numeric(data[3,4]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# MW vs W
cat("[MW vs W]")
h2_1<-as.numeric(data[2,1]); h2_1
se_1<-as.numeric(data[3,1]); se_1
h2_2<-as.numeric(data[2,5]); h2_2
se_2<-as.numeric(data[3,5]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# NE vs SE
cat("[NE vs SE]")
h2_1<-as.numeric(data[2,2]); h2_1
se_1<-as.numeric(data[3,2]); se_1
h2_2<-as.numeric(data[2,3]); h2_2
se_2<-as.numeric(data[3,3]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# NE vs SW
cat("[NE vs SW]")
h2_1<-as.numeric(data[2,2]); h2_1
se_1<-as.numeric(data[3,2]); se_1
h2_2<-as.numeric(data[2,4]); h2_2
se_2<-as.numeric(data[3,4]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# NE vs W
cat("[NE vs W]")
h2_1<-as.numeric(data[2,2]); h2_1
se_1<-as.numeric(data[3,2]); se_1
h2_2<-as.numeric(data[2,5]); h2_2
se_2<-as.numeric(data[3,5]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# SE vs SW
cat("[SE vs SW]")
h2_1<-as.numeric(data[2,3]); h2_1
se_1<-as.numeric(data[3,3]); se_1
h2_2<-as.numeric(data[2,4]); h2_2
se_2<-as.numeric(data[3,4]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# SE vs W
cat("[SE vs W]")
h2_1<-as.numeric(data[2,3]); h2_1
se_1<-as.numeric(data[3,3]); se_1
h2_2<-as.numeric(data[2,5]); h2_2
se_2<-as.numeric(data[3,5]); se_2
wald_test(h2_1, se_1, h2_2, se_2)
# SW vs W
cat("[SW vs W]")
h2_1<-as.numeric(data[2,4]); h2_1
se_1<-as.numeric(data[3,4]); se_1
h2_2<-as.numeric(data[2,5]); h2_2
se_2<-as.numeric(data[3,5]); se_2
wald_test(h2_1, se_1, h2_2, se_2)

#----------------------- Create [Heritability] Plots---------------------------#

#Create confidence intervals 
heritability_df <- heritability_df %>%
                      mutate(
                             ci_lower = h2 - 1.96 * h2_se,
                             ci_upper = h2 + 1.96 * h2_se
                            )

# set the order of X-axis
heritability_df$Sample <- factor(heritability_df$Sample, levels = c("MW", "NE", "SE", "SW", "W"))

# Draw a plot
ggplot(heritability_df, aes(x = Sample, y = h2)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
#  geom_text(aes(label = round(h2, 4)), vjust = -1.8, hjust = 0.5, size = 6) +
  geom_text(aes(y = ci_upper + 0.005, label = round(h2, 3)), size = 6) +
  labs(
    title = "Observed-Scale Heritability with 95% CI",
    subtitle = "By Region",
    x = "",
    y = "Heritability (h²)"
  ) +
  theme_minimal(base_size = 18) +
  ylim(0, max(heritability_df$ci_upper, na.rm = TRUE) + 0.01) +
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/h2_region.png"), width = 8, height = 6, dpi = 500)

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
log_files <- file.path(results_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# add birthyear range
heritability_df$byear_ranges <- byear_ranges

# View results
print(heritability_df)



#---------------------- Create [Heritability] Tables --------------------------#
## Full and Gender sample

# Transpose
heritability_birthyear <- as.data.frame(t(heritability_df))

# View results
print(heritability_birthyear)

# Table: bin4-11
latex_table_1 <- xtable(heritability_birthyear[, 1:8],
                        caption = "Heritability (Birthyear Sample: bin4–11)")

sink(file = file.path(results_dir_tables, "h2_birthyear_bin4-11.tex"))
print(latex_table_1,
      include.rownames = T,
      caption.placement = "top",
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(heritability_birthyear)))
sink()


# Table: bin12-19
latex_table_2 <- xtable(heritability_birthyear[, 9:16],
                        caption = "Heritability (Birthyear Sample: bin12–19)")

sink(file = file.path(results_dir_tables, "h2_birthyear_bin12-19.tex"))
print(latex_table_2,
      include.rownames = T,
      caption.placement = "top",
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(heritability_birthyear)))
sink()



#----------------------- Create [Heritability] Plots---------------------------#

#Create confidence intervals 
heritability_df <- heritability_df %>%
                      mutate(
                             ci_lower = h2 - 1.96 * h2_se,
                             ci_upper = h2 + 1.96 * h2_se
                            )

# set the order of X-axis
heritability_df$Sample <- factor(heritability_df$Sample, levels = paste0("bin", 4:19))

# Draw a plot
ggplot(heritability_df, aes(x = byear_ranges, y = h2)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  geom_text(aes(y = ci_upper + 0.005, label = round(h2, 3)), size = 3.3) +
  labs(
    title = "Observed-Scale Heritability with 95% CI",
    subtitle = "By Birth Year",
    x = "",
    y = "Heritability (h²)"
  ) +
  theme_minimal(base_size = 18) +
  ylim(0, max(heritability_df$ci_upper, na.rm = TRUE) + 0.01) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/h2_birthyear.png"), width = 8, height = 6, dpi = 500)

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
log_files <- file.path(results_dir, paste0(subdirs, "_h2_raw.log"))

# Apply to all files and combine results
heritability_df <- map_dfr(log_files, extract_h2_summary)

# set rowname
rownames(heritability_df) <- heritability_df$Sample

# add birthyear range
heritability_df$byear_ranges <- byear_ranges

# View results
print(heritability_df)



#---------------------- Create [Heritability] Tables --------------------------#
## Full and Gender sample

# Transpose
heritability_birthyear <- as.data.frame(t(heritability_df))

# View results
print(heritability_birthyear)

# Table: bin4-11
latex_table_1 <- xtable(heritability_birthyear,
                        caption = "Heritability (Birthyear Sample: groupingBIN1–4)")

sink(file = file.path(results_dir_tables, "h2_birthyear_groupingBINs.tex"))
print(latex_table_1,
      include.rownames = T,
      caption.placement = "top",
      booktabs = TRUE,
      hline.after = c(-1, 0, nrow(heritability_birthyear)))
sink()


n_com <- 6
# test equality of h2
wald_test <- function(h2_1, se_1, h2_2, se_2) {
  z <- (h2_1 - h2_2) / sqrt(se_1^2 + se_2^2)
  p <- 2 * pnorm(-abs(z))
  p_bon <- p*n_com
  return(c(z = z, p_value = p, p_value_bon = p_bon))
}

data <- copy(heritability_birthyear)

# 1924–1943 vs 1944–1963
cat("[1924–1943 vs 1944–1963]")
h2_1<-as.numeric(data[2,1])
 h2_1
se_1<-as.numeric(data[3,1])
 se_1
h2_2<-as.numeric(data[2,2])
 h2_2
se_2<-as.numeric(data[3,2])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)
cat("\n")
# 1924–1943 vs 1964–1983
cat("[1924–1943 vs 1964–1983]")
h2_1<-as.numeric(data[2,1])
 h2_1
se_1<-as.numeric(data[3,1])
 se_1
h2_2<-as.numeric(data[2,3])
 h2_2
se_2<-as.numeric(data[3,3])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)
cat("\n")
# 1924–1943 vs 1984–2003
cat("[1924–1943 vs 1984–2003]")
h2_1<-as.numeric(data[2,1])
 h2_1
se_1<-as.numeric(data[3,1])
 se_1
h2_2<-as.numeric(data[2,4])
 h2_2
se_2<-as.numeric(data[3,4])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)
cat("\n")
# 1944–1963 vs 1964–1983
cat("[1944–1963 vs 1964–1983]")
h2_1<-as.numeric(data[2,2])
 h2_1
se_1<-as.numeric(data[3,2])
 se_1
h2_2<-as.numeric(data[2,3])
 h2_2
se_2<-as.numeric(data[3,3])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)
cat("\n")
# 1944–1963 vs 1984–2003
cat("[1944–1963 vs 1984–2003]")
h2_1<-as.numeric(data[2,2])
 h2_1
se_1<-as.numeric(data[3,2])
 se_1
h2_2<-as.numeric(data[2,4])
 h2_2
se_2<-as.numeric(data[3,4])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)
cat("\n")
# 1964–1983 vs 1984–2003
cat("[1964–1983 vs 1984–2003]")
h2_1<-as.numeric(data[2,3])
 h2_1
se_1<-as.numeric(data[3,3])
 se_1
h2_2<-as.numeric(data[2,4])
 h2_2
se_2<-as.numeric(data[3,4])
 se_2
wald_test(h2_1, se_1, h2_2, se_2)
cat("\n")



#----------------------- Create [Heritability] Plots---------------------------#

#Create confidence intervals 
heritability_df <- heritability_df %>%
                      mutate(
                             ci_lower = h2 - 1.96 * h2_se,
                             ci_upper = h2 + 1.96 * h2_se
                            )

# set the order of X-axis
heritability_df$Sample <- factor(heritability_df$Sample, levels = paste0("groupingBIN", 1:4))

# Draw a plot
ggplot(heritability_df, aes(x = byear_ranges, y = h2)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  geom_text(aes(y = ci_upper + 0.005, label = round(h2, 3)), size = 6) +
  labs(
    title = "Observed-Scale Heritability with 95% CI",
    subtitle = "By Birth Year",
    x = "",
    y = "Heritability (h²)"
  ) +
  theme_minimal(base_size = 18) +
  ylim(0, max(heritability_df$ci_upper, na.rm = TRUE) + 0.01) +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 18, color = "gray40"),
        panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5)
       )

 # Save the plot
 ggsave(paste0(results_dir_figures,"/h2_birthyear_groupingBINs.png"), width = 8, height = 6, dpi = 500)




