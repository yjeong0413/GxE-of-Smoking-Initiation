R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

cat("Set directories\n")
start_dir="/YOUR PATH HERE/"
#sumstat_dir=paste0(start_dir,"/",ANCESTRY)

#root_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN")
#dest_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC")
annotation_dir=paste0(start_dir,"/WORKPLACE/ldsc/Annotation_Files_Clean")

ref="hg19"
reference_dir=paste0("/YOUR PATH HERE/1000G/",ref,"/clean")

#MA_dir=paste0(dest_dir,"/Meta_Analysis")


cat("\n-------------------------------------------------------------------------------\n")

cat("\nLoad Libraries\n")
cat("\n-----------------------------------------------------------------------------------\n")
cat("*Load [data.table]\n")
library(data.table)
cat("\n")
cat("*Load [dplyr]\n")
library(dplyr)
cat("\n")
cat("*Load [stringr]\n")
library(stringr)
cat("\n")
cat("*Load [ggplot2]\n")
library(ggplot2)
cat("\n")
cat("*Load [qqman]\n")
library(qqman)
cat("\n")
cat("*Load [tidyr]\n")
library(tidyr)

cat("\n-----------------------------------------------------------------------------------\n")


#ancestries = c("european", "east_asian", "african_american", "latino")
ancestry_map <-c(
                 EUR = "european",
                 EAS = "east_asian",
                 AMR = "latino",
                 AFR = "african_american"
                )

print(names(ancestry_map))
print(ancestry_map[["EUR"]])
print(ancestry_map[["EAS"]])
print(ancestry_map[["AFR"]])
print(ancestry_map[["AMR"]])

# ancestry names
ancestries <- names(ancestry_map)
ancestries

# ancestry pairs
anc_pairs <- combn(unname(ancestry_map), 2, simplify = FALSE)
anc_pairs


mods <- c(
  "female",
  "male",
  paste0("groupingBIN",1:4)
)
mods

mod_labels <- c(
  female = "Female",
  male = "Male",
  groupingBIN1 = "Silent Generation",
  groupingBIN2 = "Baby Boomers",
  groupingBIN3 = "Generation X",
  groupingBIN4 = "Millennials"
)
mod_labels

all_results <- list()

for(mod in mods){

  cat("\n")
  cat(paste0("\n*===================================== ", mod, " ================================*\n"))

  for(pair in anc_pairs){

    anc1 <- pair[1]
    anc2 <- pair[2]

    ANC <- paste0(anc1, "_", anc2)

    cat("\n")
    cat(paste0("\n*------------------------------------- ", ANC, " --------------------------------*\n"))

    file_path <- paste0(
      start_dir,
      "/WORKPLACE/ldsc/Popcorn/Output/Step2_Results/",
      ref, "/", mod, "/",
      ANC, "_ge_results_", mod, "_regression.txt"
    )

    temp <- data.frame(
      moderator = mod,
      anc1 = anc1,
      anc2 = anc2,
      file_path = file_path,
      file_exists = file.exists(file_path),

      h2_1 = NA_real_,
      h2_1_se = NA_real_,
      h2_1_p = NA_real_,

      h2_2 = NA_real_,
      h2_2_se = NA_real_,
      h2_2_p = NA_real_,

      ge = NA_real_,
      ge_se = NA_real_,
      ge_z = NA_real_,
      ge_p = NA_real_
    )

    if(!file.exists(file_path)){

      message("Missing file: ", file_path)

      all_results[[length(all_results) + 1]] <- temp
      next
    } # if

    dat <- tryCatch(
      fread(file_path),
      error = function(e){
        warning("Could not read file: ", file_path, "\n", e$message)
        return(NULL)
      } #function: error
    ) # tryCatch

    if(is.null(dat)){

      all_results[[length(all_results) + 1]] <- temp
      next
    } #if

    if(nrow(dat) == 0){

      warning("Empty file: ", file_path)

      all_results[[length(all_results) + 1]] <- temp
      next
    } #if

    names(dat)[1] <- "parameter"
    names(dat) <- gsub(" ", "_", names(dat)) # replace spaces with _
    names(dat) <- gsub("[()]", "", names(dat)) # remove () and []

    cat("\n[dat]: \n")
    cat("\nReturn first 10 rows:\n")
    print(head(dat, n = 10))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(dat))

    required_cols <- c("parameter", "Val_obs", "SE", "Z", "P_Z")

    missing_cols <- setdiff(required_cols, names(dat))

    if(length(missing_cols) > 0){

      warning(
        "Missing required columns in: ",
        file_path,
        "\nMissing columns: ",
        paste(missing_cols, collapse = ", ")
      )

      all_results[[length(all_results) + 1]] <- temp
      next
    }

    if(sum(dat$parameter == "h1^2") == 1){
      temp$h2_1 <- dat[parameter == "h1^2", Val_obs]
      temp$h2_1_se <- dat[parameter == "h1^2", SE]
      temp$h2_1_p <- dat[parameter == "h1^2", P_Z]
    } else {
      warning("Missing or duplicated h1^2 row in: ", file_path)
    }

    if(sum(dat$parameter == "h2^2") == 1){
      temp$h2_2 <- dat[parameter == "h2^2", Val_obs]
      temp$h2_2_se <- dat[parameter == "h2^2", SE]
      temp$h2_2_p <- dat[parameter == "h2^2", P_Z]
    } else {
      warning("Missing or duplicated h2^2 row in: ", file_path)
    }

    if(sum(dat$parameter == "pge") == 1){
      temp$ge <- dat[parameter == "pge", Val_obs]
      temp$ge_se <- dat[parameter == "pge", SE]
      temp$ge_z <- dat[parameter == "pge", Z]
      temp$ge_p <- dat[parameter == "pge", P_Z]
    } else {
      warning("Missing or duplicated pge row in: ", file_path)
    }

    all_results[[length(all_results) + 1]] <- temp

    cat("\n[temp]: \n")
    print(temp)

  }
}

#combine results
results <- bind_rows(all_results)
 cat("\n[results]: \n")
 cat("\nReturn first 10 rows:\n")
 print(head(results, n = 10))
 cat("\nReturn last 10 rows:\n")
 print(tail(results, n = 10))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(results))
 cat("\nNumber of missing files:\n")
 print(sum(!results$file_exists))
 cat("\nRows with missing ge:\n")
 print(results %>% filter(is.na(ge)))

all_results <- list()

for(mod in mods){
    
  cat("\n")
  cat(paste0("\n*===================================== ", mod, " ================================*\n"))

  for(pair in anc_pairs){

    anc1 <- pair[1]
    anc2 <- pair[2]

    ANC <- paste0(anc1,"_",anc2)

  cat("\n")
  cat(paste0("\n*------------------------------------- ", ANC, " --------------------------------*\n"))

    file_path <- paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Output/Step2_Results/",
                        ref,"/",mod,"/",ANC,"_ge_results_",mod,"_regression.txt")

    if(!file.exists(file_path)){
      message("Missing: ", file_path)
      next
    }

    dat <- fread(file_path)

    names(dat)[1] <- "parameter"
    names(dat) <- gsub(" ", "_", names(dat))
    names(dat) <- gsub("[()]", "", names(dat))
   
    cat("\n[dat]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(dat, n = 10))    
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(dat))

    temp <- data.frame(
      moderator = mod,
      anc1 = anc1,
      anc2 = anc2,

      h2_1 = dat[parameter == "h1^2", Val_obs],
      h2_1_se = dat[parameter == "h1^2", SE],
      h2_1_p = dat[parameter == "h1^2", P_Z],

      h2_2 = dat[parameter == "h2^2", Val_obs],
      h2_2_se = dat[parameter == "h2^2", SE],
      h2_2_p = dat[parameter == "h2^2", P_Z],

      ge = dat[parameter == "pge", Val_obs],
      ge_se = dat[parameter == "pge", SE],
      ge_z = dat[parameter == "pge", Z],
      ge_p = dat[parameter == "pge", P_Z]
    )

    all_results[[length(all_results) + 1]] <- temp

    cat("\n[all_results]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(all_results, n = 10))    
     
  } #for: pair
} #for: mod

results <- bind_rows(all_results)

    cat("\n[results]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(results, n = 10))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(results, n = 10))    
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(results))


reverse_map <- setNames(names(ancestry_map), ancestry_map)
reverse_map

panel_order <- c(
  "Female",
  "Silent Generation",
  "Baby Boomers",
  "Male",
  "Generation X",
  "Millennials"
)

plot_df <- results %>%
  mutate(
    anc1_short = reverse_map[anc1],
    anc2_short = reverse_map[anc2],
    moderator_label = mod_labels[moderator],
    moderator_label = factor(moderator_label, levels = panel_order),
    x = match(anc1_short, ancestries),
    y = match(anc2_short, ancestries),
    sig = ifelse(!is.na(ge_p) & ge_p < 0.05, "*", ""),
    ge_label = ifelse(is.na(ge), "-", paste0(sprintf("%.2f", ge), sig)),
    h2_label = ifelse(
      is.na(h2_1) | is.na(h2_2),
      "",
      paste0(anc1_short, ": h²=", sprintf("%.3f", h2_1), "\n",
             anc2_short, ": h²=", sprintf("%.3f", h2_2))
    )
  )

grid_df <- expand_grid(
  moderator_label = factor(panel_order, levels = panel_order),
  x = 1:4,
  y = 1:4
) %>%
  left_join(plot_df, by = c("moderator_label", "x", "y")) %>%
  mutate(
    cell_type = case_when(
      y > x ~ "lower",
      y == x ~ "diag",
      TRUE ~ "upper"
    ),
    fill_value = ifelse(cell_type == "lower", ge, NA_real_),
    ge_label = case_when(
      cell_type != "lower" ~ "",
      is.na(ge_label) ~ "-",
      TRUE ~ ge_label
    ),
    h2_label = case_when(
      cell_type != "lower" ~ "",
      is.na(h2_label) ~ "",
      TRUE ~ h2_label
    )
  ) %>%
  filter(cell_type == "lower")


p <- ggplot(grid_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = fill_value),
    color = "white", linewidth = 1) +
  geom_text(aes(label = ge_label),
    size = 5, nudge_y = 0.12) +
  geom_text(aes(label = h2_label),
    size = 3.5, nudge_y = -0.22) +
  scale_x_continuous(
    breaks = 1:3,
    labels = ancestries[1:3],
    position = "bottom",
    limits = c(0.5, 3.5),
    expand = c(0, 0)) +
  scale_y_reverse(
    breaks = 2:4,
    labels = ancestries[2:4],
    limits = c(4.5, 1.5),
    expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "white", # below 1
    mid = "red", # exactly 1
    high = "red", # above 1
    midpoint = 1,
    limits = c(0.4, 1.6),
    na.value = "grey95",
    name = expression(rho[ge])) +
  facet_wrap(~ moderator_label, ncol = 3) +
  coord_equal(clip = "off") +
  labs(x = NULL, y = NULL,
    title = expression(
      "Cross-ancestry genetic effect correlations (" *
        rho[ge] *") and SNP heritability (" * h^2 * ")")) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.caption = element_text(hjust = 0, size = 9),
    axis.text = element_text(color = "black"),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.spacing = unit(1.1, "lines")
  )

print(p)

# Save
file_path <- paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Output/Step2_Results/",
                    ref,"/Figures/popcorn_ge_h2_results_reg_plot.png")
ggsave(filename = file_path,
  plot = p, width = 14, height = 10, dpi = 300)


