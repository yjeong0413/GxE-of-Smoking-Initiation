R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#ANCESTRY="european"

cat("Set directories\n")
start_dir="/YOUR PATH HERE/"
root_dir=paste0(start_dir,"/WORKPLACE/ldsc/Cross-Ancestry_Meta-Analysis")
#dest_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC")

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

cat("\n-----------------------------------------------------------------------------------\n")


# Define subgroup
type <- "birthyear"
pair <- "allancestry"
#subgroups <- paste0("generation",1:4)
#subgroups

subgroups_raw <- paste0("birthyear_bin",4:7)
subgroups_raw
type2 <- "generation1"

# empty container
data_list <- list()

for (SUBGROUP in subgroups_raw) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # set file path
    raw_data_dir=paste0(root_dir,"/",type,
                        "/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned.txt")
    
    # read data
    data_temp <- fread(raw_data_dir, sep="\t", header=T)
    
    # keep only needed columns
    data <- data_temp %>% select(MarkerName, N_tot, N_cas, N_con)
    
    # add subgroup lavel 
    data[, subgroup := SUBGROUP]
    
    # store
    data_list[[SUBGROUP]] <- data
    
      cat(paste0("\nReturn first 5 rows: \n"))
      print(head(data, n = 5))    
      cat(paste0("\nReturn last 5 rows: \n"))
      print(tail(data, n = 5))
      cat("\nSize of table: (# of rows, # of columns)\n")
      print(dim(data))
    
    cat("\n===========================================================================\n")
} # for: SUBGROUP

# Stack all ancestry files vertically
data_all <- rbindlist(data_list, use.names = TRUE, fill = TRUE) # match columns by column names, not by position

# Sum sample sizes by snpid
data_merged <- data_all[, .(
                            N_tot = sum(N_tot, na.rm = F), # NA if anything is missing
                            N_cas = sum(N_cas, na.rm = F),
                            N_con = sum(N_con, na.rm = F)), 
                            by = MarkerName]

# keep final columns
data_final <- data_merged %>% select(MarkerName, N_tot, N_cas, N_con)

#assign data
data_name <- paste0("data_sample_size_", type2)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


SUBGROUP <- type2

    cat("\n------------------------------[",SUBGROUP,"]------------------------------\n")

    # Read generation1-4 meta-analysis results
    file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults1.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("")
    cat("1) Read the data.\n")
       MA_results_temp <- fread(file_path, sep = "\t", header = T)

        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_temp, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_temp, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_temp))

    # Keep SNPs that are availble across all ancestries
    # check [MA_results_temp$Direction]  
    temp <- MA_results_temp$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[BEFORE change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Keep only SNPs present in all studies
    MA_results <- MA_results_temp %>%
                            filter(!grepl("\\?",Direction))

    # Check [MA_results$Direction]  
    temp <- MA_results$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[AFTER change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Merge with [data_sample_size_XXXX]
    MA_results_cleaned <- left_join (MA_results, get(paste0("data_sample_size_",SUBGROUP)), 
                                     by = c("MarkerName"))
        cat("\n[MA_results_cleaned]\n")
        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_cleaned, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_cleaned, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_cleaned))
   
    # Save
     #Save the data
    file_path <- paste0(root_dir, "/",
                        type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt")
        cat("\nSave the results in [METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt].\n")
        write.table(MA_results_cleaned, 
                    file=file_path,
                    sep="\t", quote=F, row.names=F, col.names=T)

    cat("\n---------------------------------------------------------------------------\n")



subgroups_raw <- paste0("birthyear_bin",8:11)
subgroups_raw
type2 <- "generation2"

# empty container
data_list <- list()

for (SUBGROUP in subgroups_raw) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # set file path
    raw_data_dir=paste0(root_dir,"/",type,
                        "/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned.txt")
    
    # read data
    data_temp <- fread(raw_data_dir, sep="\t", header=T)
    
    # keep only needed columns
    data <- data_temp %>% select(MarkerName, N_tot, N_cas, N_con)
    
    # add subgroup lavel 
    data[, subgroup := SUBGROUP]
    
    # store
    data_list[[SUBGROUP]] <- data
    
      cat(paste0("\nReturn first 5 rows: \n"))
      print(head(data, n = 5))    
      cat(paste0("\nReturn last 5 rows: \n"))
      print(tail(data, n = 5))
      cat("\nSize of table: (# of rows, # of columns)\n")
      print(dim(data))
    
    cat("\n===========================================================================\n")
} # for: SUBGROUP

# Stack all ancestry files vertically
data_all <- rbindlist(data_list, use.names = TRUE, fill = TRUE) # match columns by column names, not by position

# Sum sample sizes by snpid
data_merged <- data_all[, .(
                            N_tot = sum(N_tot, na.rm = F), # NA if anything is missing
                            N_cas = sum(N_cas, na.rm = F),
                            N_con = sum(N_con, na.rm = F)), 
                            by = MarkerName]

# keep final columns
data_final <- data_merged %>% select(MarkerName, N_tot, N_cas, N_con)

#assign data
data_name <- paste0("data_sample_size_", type2)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


SUBGROUP <- type2

    cat("\n------------------------------[",SUBGROUP,"]------------------------------\n")

    # Read generation1-4 meta-analysis results
    file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults1.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("")
    cat("1) Read the data.\n")
       MA_results_temp <- fread(file_path, sep = "\t", header = T)

        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_temp, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_temp, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_temp))

    # Keep SNPs that are availble across all ancestries
    # check [MA_results_temp$Direction]  
    temp <- MA_results_temp$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[BEFORE change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Keep only SNPs present in all studies
    MA_results <- MA_results_temp %>%
                            filter(!grepl("\\?",Direction))

    # Check [MA_results$Direction]  
    temp <- MA_results$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[AFTER change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Merge with [data_sample_size_XXXX]
    MA_results_cleaned <- left_join (MA_results, get(paste0("data_sample_size_",SUBGROUP)), 
                                     by = c("MarkerName"))
        cat("\n[MA_results_cleaned]\n")
        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_cleaned, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_cleaned, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_cleaned))
   
    # Save
     #Save the data
    file_path <- paste0(root_dir, "/",
                        type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt")
        cat("\nSave the results in [METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt].\n")
        write.table(MA_results_cleaned, 
                    file=file_path,
                    sep="\t", quote=F, row.names=F, col.names=T)

    cat("\n---------------------------------------------------------------------------\n")



subgroups_raw <- paste0("birthyear_bin",12:15)
subgroups_raw
type2 <- "generation3"

# empty container
data_list <- list()

for (SUBGROUP in subgroups_raw) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # set file path
    raw_data_dir=paste0(root_dir,"/",type,
                        "/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned.txt")
    
    # read data
    data_temp <- fread(raw_data_dir, sep="\t", header=T)
    
    # keep only needed columns
    data <- data_temp %>% select(MarkerName, N_tot, N_cas, N_con)
    
    # add subgroup lavel 
    data[, subgroup := SUBGROUP]
    
    # store
    data_list[[SUBGROUP]] <- data
    
      cat(paste0("\nReturn first 5 rows: \n"))
      print(head(data, n = 5))    
      cat(paste0("\nReturn last 5 rows: \n"))
      print(tail(data, n = 5))
      cat("\nSize of table: (# of rows, # of columns)\n")
      print(dim(data))
    
    cat("\n===========================================================================\n")
} # for: SUBGROUP

# Stack all ancestry files vertically
data_all <- rbindlist(data_list, use.names = TRUE, fill = TRUE) # match columns by column names, not by position

# Sum sample sizes by snpid
data_merged <- data_all[, .(
                            N_tot = sum(N_tot, na.rm = F), # NA if anything is missing
                            N_cas = sum(N_cas, na.rm = F),
                            N_con = sum(N_con, na.rm = F)), 
                            by = MarkerName]

# keep final columns
data_final <- data_merged %>% select(MarkerName, N_tot, N_cas, N_con)

#assign data
data_name <- paste0("data_sample_size_", type2)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


SUBGROUP <- type2

    cat("\n------------------------------[",SUBGROUP,"]------------------------------\n")

    # Read generation1-4 meta-analysis results
    file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults1.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("")
    cat("1) Read the data.\n")
       MA_results_temp <- fread(file_path, sep = "\t", header = T)

        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_temp, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_temp, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_temp))

    # Keep SNPs that are availble across all ancestries
    # check [MA_results_temp$Direction]  
    temp <- MA_results_temp$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[BEFORE change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Keep only SNPs present in all studies
    MA_results <- MA_results_temp %>%
                            filter(!grepl("\\?",Direction))

    # Check [MA_results$Direction]  
    temp <- MA_results$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[AFTER change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Merge with [data_sample_size_XXXX]
    MA_results_cleaned <- left_join (MA_results, get(paste0("data_sample_size_",SUBGROUP)), 
                                     by = c("MarkerName"))
        cat("\n[MA_results_cleaned]\n")
        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_cleaned, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_cleaned, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_cleaned))
   
    # Save
     #Save the data
    file_path <- paste0(root_dir, "/",
                        type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt")
        cat("\nSave the results in [METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt].\n")
        write.table(MA_results_cleaned, 
                    file=file_path,
                    sep="\t", quote=F, row.names=F, col.names=T)

    cat("\n---------------------------------------------------------------------------\n")



subgroups_raw <- paste0("birthyear_bin",16:19)
subgroups_raw
type2 <- "generation4"

# empty container
data_list <- list()

for (SUBGROUP in subgroups_raw) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # set file path
    raw_data_dir=paste0(root_dir,"/",type,
                        "/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned.txt")
    
    # read data
    data_temp <- fread(raw_data_dir, sep="\t", header=T)
    
    # keep only needed columns
    data <- data_temp %>% select(MarkerName, N_tot, N_cas, N_con)
    
    # add subgroup lavel 
    data[, subgroup := SUBGROUP]
    
    # store
    data_list[[SUBGROUP]] <- data
    
      cat(paste0("\nReturn first 5 rows: \n"))
      print(head(data, n = 5))    
      cat(paste0("\nReturn last 5 rows: \n"))
      print(tail(data, n = 5))
      cat("\nSize of table: (# of rows, # of columns)\n")
      print(dim(data))
    
    cat("\n===========================================================================\n")
} # for: SUBGROUP

# Stack all ancestry files vertically
data_all <- rbindlist(data_list, use.names = TRUE, fill = TRUE) # match columns by column names, not by position

# Sum sample sizes by snpid
data_merged <- data_all[, .(
                            N_tot = sum(N_tot, na.rm = F), # NA if anything is missing
                            N_cas = sum(N_cas, na.rm = F),
                            N_con = sum(N_con, na.rm = F)), 
                            by = MarkerName]

# keep final columns
data_final <- data_merged %>% select(MarkerName, N_tot, N_cas, N_con)

#assign data
data_name <- paste0("data_sample_size_", type2)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


SUBGROUP <- type2

    cat("\n------------------------------[",SUBGROUP,"]------------------------------\n")

    # Read generation1-4 meta-analysis results
    file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults1.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("")
    cat("1) Read the data.\n")
       MA_results_temp <- fread(file_path, sep = "\t", header = T)

        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_temp, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_temp, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_temp))

    # Keep SNPs that are availble across all ancestries
    # check [MA_results_temp$Direction]  
    temp <- MA_results_temp$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[BEFORE change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Keep only SNPs present in all studies
    MA_results <- MA_results_temp %>%
                            filter(!grepl("\\?",Direction))

    # Check [MA_results$Direction]  
    temp <- MA_results$Direction
     ## barplot
     barplot(table(temp))
     ## summary
    cat("\n[AFTER change]: \n")
     cat("Max.: ", max(temp, na.rm = TRUE), "\n")
     cat("Min.: ", min(temp, na.rm = TRUE), "\n")
     cat("Length of non-empty: ", sum(!is.na(temp)), "\n")
     print(table(temp))

    # Merge with [data_sample_size_XXXX]
    MA_results_cleaned <- left_join (MA_results, get(paste0("data_sample_size_",SUBGROUP)), 
                                     by = c("MarkerName"))
        cat("\n[MA_results_cleaned]\n")
        cat(paste0("\nReturn first 5 rows: \n"))
        print(head(MA_results_cleaned, n = 5))    
        cat(paste0("\nReturn last 5 rows: \n"))
        print(tail(MA_results_cleaned, n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(MA_results_cleaned))
   
    # Save
     #Save the data
    file_path <- paste0(root_dir, "/",
                        type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt")
        cat("\nSave the results in [METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned_temp.txt].\n")
        write.table(MA_results_cleaned, 
                    file=file_path,
                    sep="\t", quote=F, row.names=F, col.names=T)

    cat("\n---------------------------------------------------------------------------\n")




