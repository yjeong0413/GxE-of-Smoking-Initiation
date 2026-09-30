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
annotation_dir=paste0(start_dir,"/WORKPLACE/ldsc/Annotation_Files_Clean")
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
subgroups <- paste0("generation",1:4)
subgroups

# read annotation file combining all ancestries
file_path <- paste0(annotation_dir, "/Annotation_All.txt")

cat("1) Read the data.\n")
cat(paste0("\n: [", file_path, "]\n"))
annotation_data <- fread(file_path, sep = "\t", header = T)
    
data_name <- paste0("annotation_data")
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))

# Select [snpid], [A], and [B]: This is to combine with [MA_results]
cat("\n2) Select [snpid], [A], and [B].\n")
annotation_data_selected <- annotation_data %>% select(snpid, A, B)

data_name <- paste0("annotation_data_selected") 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))

# Extract [chr], [position], and [rsid]
cat("\n3) Extract [chr], [position], and [rsid].\n")
data_name <- paste0("annotation_data")
data_info <- get(data_name) %>% select(chr, position, snpid)

data_name <- paste0("data_info") 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


for (SUBGROUP in subgroups) {
    
    cat("\n")
    cat(paste0("\n*-----------------------------------[", SUBGROUP, "]--------------------------------*\n"))
    
    file_path <- paste0(root_dir,"/",type,
                        "/METAL_cross-ancestry_", SUBGROUP,"_allancestry_MAresults_cleaned_temp.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("1) Read the data.\n")
    cleaned_data_temp <- fread(file_path, sep = "\t", header = T)
    
    cat("\n[Before change:]\n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(cleaned_data_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(cleaned_data_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(cleaned_data_temp))

    #join
    cleaned_data <- cleaned_data_temp %>%
                        inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                        inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                        mutate(across(c(A, B), tolower))
                        
    
    
    # Save results as dynamic variable
    assign(paste0("cleaned_", SUBGROUP), cleaned_data, envir = .GlobalEnv)
    
    data_name <- paste0("cleaned_", SUBGROUP)
    
    cat("\n[After change:]\n")
    cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

   cat("\n------------------------------------------------------------------------------------------\n")
   cat("\n")
} #for


for (SUBGROUP in subgroups) {
    
    cat("\n")
    cat(paste0("\n*-----------------------------------[", SUBGROUP, "]--------------------------------*\n"))
    
    #flip sign of effects to align with the raw files.

     # Read the data
     data <- get(paste0("cleaned_", SUBGROUP))

     # Change and recreate [Effect] to [Effect_new]
     data$Effect_new <- data$Effect
     data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                               -data$Effect,
                               data$Effect_new) # keep current values

     # Change and recreate [Freq1] to [Freq1_new]
     data$Freq1_new <- data$Freq1
     data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                               1-data$Freq1,
                               data$Freq1_new) # keep current values

     # Change and recreate [MinFreq] to [MinFreq_new]
     data$MinFreq_new <- data$MinFreq
     data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                               1-data$MaxFreq,
                               data$MinFreq_new) # keep current values

     # Change and recreate [MaxFreq] to [MaxFreq_new]
     data$MaxFreq_new <- data$MaxFreq
     data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                               1-data$MinFreq,
                               data$MaxFreq_new) # keep current values

     # Select necessary variables
     data_new <- data %>%
                    select(chr, position, MarkerName, A, B, 
                           Effect_new, StdErr, "P-value", Freq1_new,
                           N_tot, N_cas, N_con)

        # Assign results as dynamic variable
        assign(paste0("cleaned_",SUBGROUP,"_new"), data_new, envir = .GlobalEnv)
        data_name <- paste0("cleaned_",SUBGROUP,"_new")

        cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
        print(head(get(data_name), n = 5))    
        cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
        print(tail(get(data_name), n = 5))
        cat("\nSize of table: (# of rows, # of columns)\n")
        print(dim(get(data_name)))

         #Save the data: _cleaned
         file_path <- paste0(root_dir, "/",
                            type,"/METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned.txt")
            cat("\nSave the results in [METAL_cross-ancestry_",SUBGROUP,"_",pair,"_MAresults_cleaned.txt].\n")
            write.table(get(data_name), 
                        file=file_path,
                        sep="\t", quote=F, row.names=F, col.names=T)

        #Save the data for FUMA: hg_38
        file_path <- paste0(root_dir, "/",
                            type,"/FUMA/Results_",type,"_",SUBGROUP,"_hg38.txt")
        cat("Save the results in [Results_",type,"_",SUBGROUP,"_hg38.txt].\n")
        write.table(get(data_name), 
                    file=file_path,
                    sep="\t", quote=F, row.names=F, col.names=T)

   cat("\n------------------------------------------------------------------------------------------\n")
   cat("\n")
} #for


type2="generation1234"
subgroups2 <- c("generation1","generation2","generation3","generation4")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")


type2="generation12"
subgroups2 <- c("generation1","generation2")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")


type2="generation13"
subgroups2 <- c("generation1","generation3")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")


type2="generation14"
subgroups2 <- c("generation1","generation4")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")


type2="generation23"
subgroups2 <- c("generation2","generation3")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")


type2="generation24"
subgroups2 <- c("generation2","generation4")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")


type2="generation34"
subgroups2 <- c("generation3","generation4")
subgroups2

# Read MA results
file_path <- paste0(root_dir, "/",type,"/METAL_cross-ancestryMA_",type2,"_MAresults1.txt")
cat(paste0("\n1. [", file_path, "]\n"))

cat("1) Read the data.\n")
   MA_results_temp <- fread(file_path, sep = "\t", header = T)


    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(MA_results_temp, n = 5))    
    cat(paste0("\nReturn last 5 rows: \n"))
    print(tail(MA_results_temp, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(MA_results_temp))

# empty container
data_list <- list()

for (SUBGROUP in subgroups2) {
    cat("\n==============================[",SUBGROUP,"]==============================\n")

    # read data
    data_temp <- get(paste0("cleaned_",SUBGROUP,"_new"))
    
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
data_name <- paste0("data_sample_size_", type)
assign(data_name, data_final)
                          
 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))
                          
cat("\n===========================================================================\n")


#1. Keep all subgroups have non-missing values
temp <- MA_results_temp$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[BEFORE change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)

# Keep only SNPs present in all studies
MA_results_temp2 <- MA_results_temp %>%
                        filter(!grepl("\\?",Direction))

# check [MA_results$Direction]  
temp <- MA_results_temp2$Direction
 ## barplot
 barplot(table(temp))
 ## summary
cat("\n[AFTER change]: \n")
 cat("Max.: ", max(temp, na.rm = TRUE), "\n")
 cat("Min.: ", min(temp, na.rm = TRUE), "\n")
 cat("Length of non-empty: ", sum(!is.na(temp)))
 table(temp)


#2. Combine with [data_sample_size_XXXX], [annotation_data_selected], [data_info]
MA_results <- MA_results_temp2 %>%
                    inner_join(get(paste0("data_sample_size_", type)), by = "MarkerName") %>%
                    inner_join(annotation_data_selected, by=c("MarkerName"="snpid")) %>%
                    inner_join(data_info, by=c("MarkerName"="snpid")) %>%
                    mutate(across(c(A, B), tolower))


 cat(paste0("\nReturn first 5 rows: \n"))
 print(head(MA_results, n = 5))    
 cat(paste0("\nReturn last 5 rows: \n"))
 print(tail(MA_results, n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(MA_results))


#3. Allele frequency check# check [MA_results_temp$Direction]  

 # Read the data
 data <- MA_results

 # Change and recreate [Effect] to [Effect_new]
 data$Effect_new <- data$Effect
 data$Effect_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           -data$Effect,
                           data$Effect_new) # keep current values

 # Change and recreate [Freq1] to [Freq1_new]
 data$Freq1_new <- data$Freq1
 data$Freq1_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$Freq1,
                           data$Freq1_new) # keep current values

 # Change and recreate [MinFreq] to [MinFreq_new]
 data$MinFreq_new <- data$MinFreq
 data$MinFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MaxFreq,
                           data$MinFreq_new) # keep current values

 # Change and recreate [MaxFreq] to [MaxFreq_new]
 data$MaxFreq_new <- data$MaxFreq
 data$MaxFreq_new <- ifelse(as.character(data$A)==as.character(data$Allele1),
                           1-data$MinFreq,
                           data$MaxFreq_new) # keep current values

# Assign results as dynamic variable
assign(paste0("MA_results_new"), data, envir = .GlobalEnv)
data_name <- paste0("MA_results_new")
 
 cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
 print(head(get(data_name), n = 5))    
 cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
 print(tail(get(data_name), n = 5))
 cat("\nSize of table: (# of rows, # of columns)\n")
 print(dim(get(data_name)))


# Save the data: _cleaned
file_path <- paste0(root_dir, "/",
                    type,"/METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt")
cat("\nSave the results in [METAL_cross-ancestryMA_",type2,"_MAresults_cleaned.txt].\n")
write.table(get(data_name), 
            file=file_path,
            sep="\t", quote=F, row.names=F, col.names=T)


# Extract/save data to use for FUMA
MA_data_new <- MA_results_new %>%
                   dplyr::select(chr, position, MarkerName, A, B, 
                                 Effect_new, StdErr, "P-value", Freq1_new,
                                 HetPVal,
                                 N_tot, N_cas, N_con)

 cat(paste0("\nReturn first 5 rows of [MA_data_new]\n"))
 print(head(MA_data_new, n = 5))    
 cat(paste0("\nReturn last 5 rows of [MA_data_new]\n"))
 print(tail(MA_data_new, n = 5))
 cat("\nSize of [MA_data_new]: (# of rows, # of columns)\n")
 print(dim(MA_data_new))
 cat("\n[chr]")
 table(MA_data_new$chr)

 #Save the data: _hg38
 file_path <- paste0(root_dir, "/",
                     type,"/FUMA/MA_results_",type2,"_hg38.txt")
 cat("\nSave the results in [MA_results_",type2,"_hg38.txt].\n")
 write.table(MA_data_new, 
             file=file_path,
             sep="\t", quote=F, row.names=F, col.names=T)



# Create 
data_MA_sigHetPVal <- MA_results_new %>%
                        filter(-log10(HetPVal) > -log10(5*10^{-8}))
    cat(paste0("\nReturn first 5 rows of [data_MA_sigHetPVal]\n"))
    print(head(data_MA_sigHetPVal, n = 5))    
    cat(paste0("\nReturn last 5 rows of [data_MA_sigHetPVal]\n"))
    print(tail(data_MA_sigHetPVal, n = 5))
    cat("\nSize of [MA_results_new]: (# of rows, # of columns)\n")
    print(dim(MA_results_new))
    cat("\nSize of [data_MA_sigHetPVal]: (# of rows, # of columns)\n")
    print(dim(data_MA_sigHetPVal))
    cat("\n")

cat("\n [data_MA_sigHetPVal$chr]\n")
    print(table(data_MA_sigHetPVal$chr))
cat("\n")



