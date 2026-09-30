R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#Define directories for "EAST ASIAN"
cat("Set ancestries\n")
ANCESTRY="east_asian"

cat("Set directories\n")
start_dir="/YOUR PATH HERE/"
snp_dir=paste0(start_dir,"/v10.2_annotation_files/",ANCESTRY,"_ac45")
sumstat_dir=paste0(start_dir,"/",ANCESTRY)
root_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN")
dest_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC")

cat("\n-------------------------------------------------------------------------------\n")


cat("\nLoad Libraries\n")
cat("\n-----------------------------------------------------------------------------------\n")
cat("*Load [data.table]\n")
library(data.table)
cat("\n")
cat("*Load [dplyr]\n")
library(dplyr)
cat("\n")
cat("*Load [qqman]\n")
library(qqman)
cat("\n-----------------------------------------------------------------------------------\n")


#Set directory

#Current directory
#cat("Current directory\n")
#setwd(paste0(destination_dir,"/ever_tobacco_user/"))
#cat(getwd())
#cat("\n-------------------------------------------------------------------------------\n")


#Function to READ summary statistics data
read_sumstats_subgroups <- function(sub, dir) {
  for (SUBGROUP in sub) {
    cat("\n")
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))
    cat(paste0("\n*---------------Start Reading: [ever_tobacco_user_", SUBGROUP, "_passed.txt]----------------*\n"))
    
    file_path <- paste0(dir, "/ever_tobacco_user_", SUBGROUP,"/ever_tobacco_user_", SUBGROUP, "_passed.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("1) Read the data.\n")
    
    sumstats <- fread(file_path, sep = "\t", skip = 1, header = FALSE)
    
    cat("- Assign column names.\n")
    colnames(sumstats) <- c("all.data.id", "src", "log.p", "effect", "stderr", "pass",
                            "im.num.0", "dose.b.0", "im.num.1", "dose.b.1",
                            "AA.0", "AB.0", "BB.0", "AA.1", "AB.1", "BB.1", 
                            "N_tot", "N_cas", "N_con")
   
     
#Create Allele Frequencies separately for case (1) and control (0)
    # Genotyped SNPs ([src]==G)
    ## Total number of A alleles (N_A)
    ### control
    sumstats$N_A0 <- NA  # initialize
    sumstats$N_A0 <- ifelse(
                            sumstats$src == "G",
                            (2*sumstats$AA.0+sumstats$AB.0),
                            sumstats$N_A0 # keep existing values.
                            ) #ifelse
    ### case
    sumstats$N_A1 <- NA  # initialize
    sumstats$N_A1 <- ifelse(
                            sumstats$src == "G",
                            (2*sumstats$AA.1+sumstats$AB.1),
                            sumstats$N_A1 # keep existing values.
                            ) #ifelse
      
    ## Total number of B alleles (N_B)
    ### control
    sumstats$N_B0 <- NA  # initialize
    sumstats$N_B0 <- ifelse(
                            sumstats$src == "G",
                            (2*sumstats$BB.0+sumstats$AB.0),
                            sumstats$N_B0 # keep existing values.
                            ) #ifelse
    ### case
    sumstats$N_B1 <- NA  # initialize
    sumstats$N_B1 <- ifelse(
                            sumstats$src == "G",
                            (2*sumstats$BB.1+sumstats$AB.1),
                            sumstats$N_B1 # keep existing values.
                            ) #ifelse

    ## Total number of alleles in controls (N_0=N_A0+N_B0) and cases (N_1=N_A1+N_B1)
    ### control
    sumstats$N_0 <- NA  # initialize
    sumstats$N_0 <- ifelse(
                            sumstats$src == "G",
                            (sumstats$N_A0+sumstats$N_B0),
                            sumstats$N_0 # keep existing values.
                            ) #ifelse
    ### case
    sumstats$N_1 <- NA  # initialize
    sumstats$N_1 <- ifelse(
                            sumstats$src == "G",
                            (sumstats$N_A1+sumstats$N_B1),
                             sumstats$N_1 # keep existing values.
                            ) #ifelse

    ## frequency of A allele ([freq.a0] and [freq.a1])
    ### control
    sumstats$freq.a0 <- NA  # initialize
    sumstats$freq.a0 <- ifelse(
                                sumstats$src == "G",
                                (sumstats$N_A0/sumstats$N_0),
                                sumstats$freq.a0 # keep existing values.
                                ) #ifelse
    ### case
    sumstats$freq.a1 <- NA  # initialize
    sumstats$freq.a1 <- ifelse(
                                sumstats$src == "G",
                                (sumstats$N_A1/sumstats$N_1),
                                sumstats$freq.a1 # keep existing values.
                                ) #ifelse

    ## frequency of B allele ([freq.b0] and [freq.b1])
    ### control
    sumstats$freq.b0 <- NA  # initialize
    sumstats$freq.b0 <- ifelse(
                                sumstats$src == "G",
                                (sumstats$N_B0/sumstats$N_0),
                                sumstats$freq.b0 # keep existing values.
                                ) #ifelse
    ### case
    sumstats$freq.b1 <- NA  # initialize
    sumstats$freq.b1 <- ifelse(
                                sumstats$src == "G",
                                (sumstats$N_B1/sumstats$N_1),
                                sumstats$freq.b1 # keep existing values.
                                ) #ifelse

    # Imputed SNPs ([src]==I)
    ## frequency of B allele (freq.b0 and freq.b1)
    ### control
    sumstats$freq.b0 <- ifelse(
                                sumstats$src == "I",
                                sumstats$dose.b.0,
                                sumstats$freq.b0 # keep existing values.
                                ) #ifelse
    ### case
    sumstats$freq.b1 <- ifelse(
                                sumstats$src == "I",
                                sumstats$dose.b.1,
                                sumstats$freq.b1 # keep existing values.
                                ) #ifelse

    ## frequency of A allele (freq.a)
    ### control
    sumstats$freq.a0 <- ifelse(
                                sumstats$src == "I",
                                (1-sumstats$freq.b0),
                                sumstats$freq.a0 # keep existing values.
                                ) #ifelse
    ### case
    sumstats$freq.a1 <- ifelse(
                                sumstats$src == "I",
                                (1-sumstats$freq.b1),
                                sumstats$freq.a1 # keep existing values.
                                ) #ifelse
      

    #Generate[MAF]
    # Genotyped SNPs ([src]==G)
    ## Total number of A alleles (N_A)
    sumstats$N_A <- NA  # initialize
    sumstats$N_A <- ifelse(
                            sumstats$src == "G",
                            ( 2*sumstats$AA.0+sumstats$AB.0
                             +2*sumstats$AA.1+sumstats$AB.1),
                             sumstats$N_A # keep existing values.
                            ) #ifelse

    ## Total number of B alleles (N_B)
    sumstats$N_B <- NA  # initialize
    sumstats$N_B <- ifelse(
                            sumstats$src == "G",
                            ( 2*sumstats$BB.0+sumstats$AB.0
                             +2*sumstats$BB.1+sumstats$AB.1),
                             sumstats$N_B # keep existing values.
                            ) #ifelse

    ## Total number of alleles (N_AB=N_A+N_B)
    sumstats$N_AB <- NA  # initialize
    sumstats$N_AB <- ifelse(
                             sumstats$src == "G",
                             (sumstats$N_A+sumstats$N_B),
                              sumstats$N_AB # keep existing values.
                             ) #ifelse

    ## frequency of A allele (freq.a)
    sumstats$freq.a <- NA  # initialize
    sumstats$freq.a <- ifelse(
                               sumstats$src == "G",
                               (sumstats$N_A/sumstats$N_AB),
                               sumstats$freq.a # keep existing values.
                              ) #ifelse

    ## frequency of B allele (freq.b)
    sumstats$freq.b <- NA  # initialize
    sumstats$freq.b <- ifelse(
                               sumstats$src == "G",
                               (sumstats$N_B/sumstats$N_AB),
                               sumstats$freq.b # keep existing values.
                               ) #ifelse
 
     # Imputed SNPs ([src]==I)
    ## frequency of B allele (freq.b)
    sumstats$freq.b <- ifelse(
                               sumstats$src == "I",
                               ((sumstats$im.num.0*sumstats$dose.b.0+sumstats$im.num.1*sumstats$dose.b.1)/(sumstats$im.num.0+sumstats$im.num.1)),
                                sumstats$freq.b # keep existing values.
                               ) #ifelse

    ## frequency of A allele (freq.a)
    sumstats$freq.a <- ifelse(
                               sumstats$src == "I",
                               (1-sumstats$freq.b),
                               sumstats$freq.a # keep existing values.
                               ) #ifelse
      
    # check whether [freq.a] and [freq.b] are all filled in.
    cat("\nCheck whether [freq.a] and [freq.b] are all filled in.\n")
    cat("-[freq.a] :",sum(is.na(sumstats$freq.a)), "\n")
    cat("-[freq.b] :",sum(is.na(sumstats$freq.b)), "\n")  
    
    # [freq.a]+[freq.b] 
    data.a <- sumstats$freq.a
    data.b <- sumstats$freq.b
    data <- data.a + data.b
    ## histogram
    hist(data)
    ## summary
    cat("Max. of [freq.a]+[freq.b]: ", max(data, na.rm = TRUE), "\n")
    cat("Min. of [freq.a]+[freq.b]: ", min(data, na.rm = TRUE), "\n")
    cat("Mean. of [freq.a]+[freq.b]: ", mean(data, na.rm = TRUE), "\n")
    cat("Median. of [freq.a]+[freq.b]: ", median(data, na.rm = TRUE), "\n")
    cat("Length of non-empty [freq.a]+[freq.b]: ", sum(!is.na(data)))
    
    #Generate [MAF]
    sumstats$MAF <- pmin(sumstats$freq.a, sumstats$freq.b)
  
    cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
    print(head(sumstats, n=5))
    cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
    print(tail(sumstats, n=5))
    cat("\nSize of a table: (# of rows, # of columns)\n")
    print(dim(sumstats))
    cat("\nCount if there are any missing values in 'all.data.id:\n'")
    print(sum(is.na(sumstats$all.data.id))) 
      
      
      
    # Save results as dynamic variable
    assign(paste0("sumstats_", SUBGROUP), sumstats, envir = .GlobalEnv)
    
    data_name <- paste0("sumstats_", SUBGROUP)
    
    cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
    print(head(get(data_name), n = 5))
    
    cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
    print(tail(get(data_name), n = 5))
    
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(sumstats))
    
    cat(paste0("\n*-------------------End Reading: [ever_tobacco_user_", SUBGROUP, "_passed.txt]--------------*\n"))
    cat("\n")
  } #for
} #function

#Restructure data
restructure_sumstats_subgroups <- function(sub, dir) {
  for (SUBGROUP in sub) {
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))
    data_name <- paste0("sumstats_", SUBGROUP)
    data_name_new <- paste0("sumstats_", SUBGROUP, "_subset")

    # keep certain columns in [sumstat_SUBGROUP]
    cat("\nkeep certain columns in [sumstat_SUBGROUP].\n")
    assign(data_name_new, subset(get(data_name), select=c("all.data.id", "src", "log.p", "effect", "stderr", 
                                                          "N_tot", "N_cas", "N_con", "pass",
                                                          "freq.a0", "freq.a1", "freq.b0", "freq.b1",
                                                          "freq.a", "freq.b", "MAF")), envir = .GlobalEnv)
    print(head(get(data_name_new),n=5))
    print(dim(get(data_name_new)))
    # remove original data to save space  
    rm(list = data_name, envir = .GlobalEnv)

    # merge [sumstat_SUBGROUP_subset] and [snp_all_merged_subset] using "all.data.id"
    cat("\nmerge [sumstat_SUBGROUP_subset] and [snp_all_merged_subset] using all.data.id.\n")
    merged <- merge(get(data_name_new), snp_all_merged_subset, by="all.data.id") # I do not specify "all=T" as I only need the matching SNP.
    print(head(merged, n=5))
    print(tail(merged, n=5))
    print(dim(merged))

    #Remove 'chr' from 'scaffold'.
    cat("\nRemove 'chr' from 'scaffold'.\n")
    merged_chr_CLEAN <- merged %>%
                            mutate(chr=substring(scaffold, 4)) # Substring start from the 4th character.
    print(head(merged_chr_CLEAN, n=5))
    print(tail(merged_chr_CLEAN, n=5))
    print(dim(merged_chr_CLEAN))

    #Separate 'Alleles' into 'Allele 1' and 'Allele 2'
    cat("\nSeparate 'Alleles' column into 'A' and 'B'.\n")
    merged_chr_allele_CLEAN <- merged_chr_CLEAN %>%
                               mutate(
                                      A=substr(alleles,1,1), #first character
                                      B=substr(alleles,nchar(alleles),nchar(alleles)) #last character
                                      ) #mutate
    print(head(merged_chr_allele_CLEAN, n=5))
    print(tail(merged_chr_allele_CLEAN, n=5))
    print(dim(merged_chr_allele_CLEAN))

    #Create [log10p].
    cat("\nCreate [log10p].\n")
    merged_chr_allele_CLEAN$log10p <- merged_chr_allele_CLEAN$log.p/log(10)

    print(head(merged_chr_allele_CLEAN, n=5))
    print(tail(merged_chr_allele_CLEAN, n=5))
    print(dim(merged_chr_allele_CLEAN))

    #Take subset of [merged_chr_allele_CLEAN]
    merged_CLEAN1 <- merged_chr_allele_CLEAN %>%
                        subset(select = c("all.data.id", "src", "chr", "position", "assay.name.all", "A", "B", 
                                          "log.p", "log10p", "effect", "stderr", 
                                          "N_tot", "N_cas", "N_con", 
                                          "freq.a0", "freq.a1", "freq.b0", "freq.b1",
                                          "freq.a", "freq.b", "MAF",
                                          "im.freq.a", "im.freq.b", "avg.rsqr", "min.rsqr", 
                                          "gt.rate", "gt.freq.a", "gt.freq.b"))
    #Assign new names
    colnames(merged_CLEAN1)<-c("all.data.id", "src", "chr", "position", "snpid", "A", "B", 
                              "log.p", "log10p", "effect", "stderr", 
                              "N_tot", "N_cas", "N_con", 
                              "freq.a0", "freq.a1", "freq.b0", "freq.b1",
                              "freq.a", "freq.b", "MAF",
                              "im.freq.a", "im.freq.b", "avg.rsqr", "min.rsqr", 
                              "gt.rate", "gt.freq.a", "gt.freq.b")


    print(head(merged_CLEAN1, n=5))
    print(tail(merged_CLEAN1, n=5))
    print(dim(merged_CLEAN1))


   ##-----------------------------------
   cat("\nCreate [MAF_old].\n")
    #Initialize [MAF_old]
    merged_CLEAN1$MAF_old <- NA

    cat("\n-Genotyped case ([src == G]).\n")
    #Genotyped case ([src == G])
    ##Replace [MAF_old] - genotyped ([src == G])
    merged_CLEAN1$MAF_old <- ifelse(
                                merged_CLEAN1$src == "G",
                                pmin(merged_CLEAN1$gt.freq.a, merged_CLEAN1$gt.freq.b),
                                merged_CLEAN1$MAF_old  # keep existing values.                            
                                ) #ifelse

    print(head(merged_CLEAN1, n=5))
    print(tail(merged_CLEAN1, n=5))
    print(dim(merged_CLEAN1))

    # [MAF_old] <- when only genotyped case is considered
    ## histogram
    hist(merged_CLEAN1$MAF_old)
    ## summary
    cat("Max. of [MAF_old]: ", max(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Min. of [MAF_old]: ", min(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Mean. of [MAF_old]: ", mean(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Median. of [MAF_old]: ", median(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Length of non-empty [MAF_old]: ", sum(!is.na(merged_CLEAN1$MAF_old)))

    cat("\n-Imputed case ([src == G]).\n")
    #Imputed case ([src == I])
    ##Replace [MAF_old] - Imputed ([src == I])
    merged_CLEAN1$MAF_old <- ifelse(
                                merged_CLEAN1$src == "I",
                                pmin(merged_CLEAN1$im.freq.a, merged_CLEAN1$im.freq.b),
                                merged_CLEAN1$MAF_old  # keep existing values (e.g., those from src == "G")                            
                                ) #ifelse

    print(head(merged_CLEAN1, n=5))
    print(tail(merged_CLEAN1, n=5))
    print(dim(merged_CLEAN1))


    # [MAF_old] <- when both genotyped and imputed cases are considered
    ## histogram
    hist(merged_CLEAN1$MAF_old)
    ## summary
    cat("Max. of [MAF_old]: ", max(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Min. of [MAF_old]: ", min(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Mean. of [MAF_old]: ", mean(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Median. of [MAF_old]: ", median(merged_CLEAN1$MAF_old, na.rm = TRUE), "\n")
    cat("Length of non-empty [MAF_old]: ", sum(!is.na(merged_CLEAN1$MAF_old)))

    #Save the data BEFORE removing SNPs by [MAF]
    cat("\nSave the results in [merged_CLEAN1.txt].\n")
    #write.table(merged_CLEAN1, file=paste0(dir,"/ever_tobacco_user_",SUBGROUP,"/merged_CLEAN1.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    cat("\n")
  ##-----------------------------------


    #Remove SNPs with MAF < 0.01
    cat("\n-Remove SNPs with MAF < 0.01.\n")
    merged_CLEAN2 <- merged_CLEAN1 %>%
                        filter(MAF>=0.01) #Automatically drops rows with MAF = NA as well (since NA >= 0.001 is NA, and filter() drops NA by default)

    print(head(merged_CLEAN2, n=5))
    print(tail(merged_CLEAN2, n=5))
    print(dim(merged_CLEAN2))

    # [MAF] <- when both genotyped and imputed cases are considered
    temp <- merged_CLEAN2$MAF
    ## histogram
    hist(temp)
    ## summary
    cat("Max. of [MAF]: ", max(temp, na.rm = TRUE), "\n")
    cat("Min. of [MAF]: ", min(temp, na.rm = TRUE), "\n")
    cat("Mean. of [MAF]: ", mean(temp, na.rm = TRUE), "\n")
    cat("Median. of [MAF]: ", median(temp, na.rm = TRUE), "\n")
    cat("Length of non-empty [MAF]: ", sum(!is.na(temp)))


    #Take subset of [merged_CLEAN]
    merged_CLEAN2 <- merged_CLEAN2 %>%
                        subset(select = c("all.data.id","src",
                                          "chr", "position", "snpid", "A", "B", 
                                          "log.p", "log10p", "effect", "stderr", 
                                          "N_tot", "N_cas", "N_con",
                                          "freq.a0", "freq.a1", "freq.b0", "freq.b1",
                                          "freq.a", "freq.b", "MAF", "MAF_old"
                                          ) # select
                              ) # subset


    #Save the data AFTER removing SNPs by [MAF]
    cat("\nSave the results in [merged_CLEAN2.txt].\n")
    write.table(merged_CLEAN2, file=paste0(dir,"/ever_tobacco_user_",SUBGROUP,"/merged_CLEAN2.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    cat("\n")

    #Export the data in the global environment (?)
    data_name_export <- paste0("merged_CLEAN2_", SUBGROUP)
    assign(data_name_export, merged_CLEAN2, envir = .GlobalEnv)

    cat(paste0("\n*------------------------------------------------------------------------------------*\n"))
     
  } #for
} #function

# Generate data for munge
munge_sumstats_subgroups <- function(sub, dir) {
  for (SUBGROUP in sub) {
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))

    #read the data
    merged_CLEAN2 <- get(paste0("merged_CLEAN2_", SUBGROUP))


    #Keep 'chr1-22' only.
    cat("\nKeep 'chr1-22' only.\n")
    merged_CLEAN3_temp <- merged_CLEAN2 %>%
                        filter(chr %in% as.character(1:22)) # Keep chr1-22 only.

    print(head(merged_CLEAN3_temp, n=5))
    print(tail(merged_CLEAN3_temp, n=5))
    print(dim(merged_CLEAN3_temp))


    #Create raw p-value.
    cat("\nCreate raw p-value.\n")
    merged_CLEAN3_temp$pval <- exp(merged_CLEAN3_temp$log.p)

    print(head(merged_CLEAN3_temp, n=5))
    print(tail(merged_CLEAN3_temp, n=5))
    print(dim(merged_CLEAN3_temp))

    # [pval] <- raw p-value
    temp <- merged_CLEAN3_temp$pval
    ## histogram
    hist(temp)
    ## summary
    cat("Max. of [pval]: ", max(temp, na.rm = TRUE), "\n")
    cat("Min. of [pval]: ", min(temp, na.rm = TRUE), "\n")
    cat("Mean. of [pval]: ", mean(temp, na.rm = TRUE), "\n")
    cat("Median. of [pval]: ", median(temp, na.rm = TRUE), "\n")
    cat("Length of non-empty [pval]: ", sum(!is.na(temp)))

     # [log.p] <- natural log p
    temp <- merged_CLEAN3_temp$log.p
    ## histogram
    hist(temp)
    ## summary
    cat("Max. of [log.p]: ", max(temp, na.rm = TRUE), "\n")
    cat("Min. of [log.p]: ", min(temp, na.rm = TRUE), "\n")
    cat("Mean. of [log.p]: ", mean(temp, na.rm = TRUE), "\n")
    cat("Median. of [log.p]: ", median(temp, na.rm = TRUE), "\n")
    cat("Length of non-empty [log.p]: ", sum(!is.na(temp)))

     
    #Take subset of [merged_CLEAN3_temp].
    cat("\nTake subset of [merged_CLEAN3_temp].\n")
    merged_CLEAN3 <- merged_CLEAN3_temp %>%
                        subset(select = c("chr", "position", "snpid", "A", "B", 
                                          "log.p", "pval", "effect", "stderr", 
                                          "N_tot", "N_cas", "N_con",
                                          "MAF") #select
                              ) #subset

    print(head(merged_CLEAN3, n=5))
    print(tail(merged_CLEAN3, n=5))
    print(dim(merged_CLEAN3))

    #Save the data
    cat("Save the results in [merged_CLEAN3.txt].\n")
    write.table(merged_CLEAN3, file=paste0(dir,"/ever_tobacco_user_",SUBGROUP,"/merged_CLEAN3.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    cat("\n")
    cat(paste0("\n*------------------------------------------------------------------------------------*\n"))
      
 } #for
} #function


#Generate data including [allele frequencies] from [merged_CLEAN2.txt] ([merged_CLEAN4.txt])
metaanalysis_sumstats_subgroups <- function(sub, dir) {
  for (SUBGROUP in sub) {
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))

    #read the data
    merged_CLEAN2 <- get(paste0("merged_CLEAN2_", SUBGROUP))

    # Create [merged_CLEAN4_temp]
    merged_CLEAN4_temp <- merged_CLEAN2
      
    print(head(merged_CLEAN4_temp, n=5))
    print(tail(merged_CLEAN4_temp, n=5))
    print(dim(merged_CLEAN4_temp))

     #[src]
     print(table(merged_CLEAN4_temp$src))
      
     #[chr]
     print(table(merged_CLEAN4_temp$chr))
      
    #Create raw p-value.
    cat("\nCreate raw p-value.\n")
    merged_CLEAN4_temp$pval <- exp(merged_CLEAN4_temp$log.p)

    print(head(merged_CLEAN4_temp, n=5))
    print(tail(merged_CLEAN4_temp, n=5))
    print(dim(merged_CLEAN4_temp))

    # [pval] <- raw p-value
    temp <- merged_CLEAN4_temp$pval
    ## histogram
    hist(temp)
    ## summary
    cat("Max. of [pval]: ", max(temp, na.rm = TRUE), "\n")
    cat("Min. of [pval]: ", min(temp, na.rm = TRUE), "\n")
    cat("Mean. of [pval]: ", mean(temp, na.rm = TRUE), "\n")
    cat("Median. of [pval]: ", median(temp, na.rm = TRUE), "\n")
    cat("Length of non-empty [pval]: ", sum(!is.na(temp)))

     # [log.p] <- natural log p
    temp <- merged_CLEAN4_temp$log.p
    ## histogram
    hist(temp)
    ## summary
    cat("Max. of [log.p]: ", max(temp, na.rm = TRUE), "\n")
    cat("Min. of [log.p]: ", min(temp, na.rm = TRUE), "\n")
    cat("Mean. of [log.p]: ", mean(temp, na.rm = TRUE), "\n")
    cat("Median. of [log.p]: ", median(temp, na.rm = TRUE), "\n")
    cat("Length of non-empty [log.p]: ", sum(!is.na(temp)))

     
    #Take subset of [merged_CLEAN3_temp].
    cat("\nTake subset of [merged_CLEAN3_temp].\n")
    merged_CLEAN4 <- merged_CLEAN4_temp %>%
                        subset(select = c("chr", "position", "snpid", "A", "B", 
                                          "log.p", "pval", "effect", "stderr", 
                                          "N_tot", "N_cas", "N_con",
                                          "freq.a0", "freq.a1", "freq.b0", "freq.b1",
                                          "freq.a", "freq.b",
                                          "MAF") #select
                              ) #subset

    print(head(merged_CLEAN4, n=5))
    print(tail(merged_CLEAN4, n=5))
    print(dim(merged_CLEAN4))

    #Save the data
    cat("Save the results in [merged_CLEAN4.txt].\n")
    write.table(merged_CLEAN4, file=paste0(dir,"/ever_tobacco_user_",SUBGROUP,"/merged_CLEAN4.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    cat("\n")
    cat(paste0("\n*------------------------------------------------------------------------------------*\n"))
      
 } #for
} #function


#Read txt files
cat("\n*---------------------------------- Merged SNP file -------------------------------------*\n")
cat("\n------------------------------------------------------------------------------------------\n")

cat("[snp_all_merged.txt]\n")                                                                 
cat("\n1) Read the data.\n")                                                                    
snp_all_merged <- fread(paste0(root_dir,"/snp/Merged/snp_all_merged.txt"), sep="\t", skip=1, header=F)
cat("\n-Assign column names.\n")
colnames(snp_all_merged)<-c("gt.data.id","im.data.id","all.data.id",
                            "assay.name.all","scaffold","position","alleles",
                            "ploidy", "cytoband", "gene.context", 
                            "is.v1", "is.v2", "is.v3", "is.v4", "is.v5", 
                            "h550", "omni", "strand",
                            "assay.name.im", "im.freq.a", "im.freq.b", "avg.rsqr", "min.rsqr", 
                            "p.batch", "qc.mask",
                            "assay.name.gt", "gt.rate", "gt.freq.a", "gt.freq.b", "hw.p.value", "p.date"
                           )                                              
cat("\nReturn first 5 rows of [snp_all_merged.txt]\n")
print(head(snp_all_merged, n=5))
cat("\nReturn last 5 rows of [snp_all_merged.txt]\n")
print(tail(snp_all_merged, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")                                         
print(dim(snp_all_merged))
cat("\nCount if there are any missing values in 'gt.data.id':\n")                             
print(sum(is.na(snp_all_merged$all.data.id)))

# keep certain columns in [snp_all_merged]
snp_all_merged_subset <- subset(snp_all_merged, select=c("all.data.id", "assay.name.all","scaffold","position","alleles", "im.freq.a", "im.freq.b", "avg.rsqr", "min.rsqr", "gt.rate", "gt.freq.a", "gt.freq.b"))
print(head(snp_all_merged_subset,n=5))
print(dim(snp_all_merged_subset))

#Remove data to save memory
rm(snp_all_merged)

#Define subgroup
subgroups <- paste0("birthyear_bin", 4:8)

#Read summary statistics data
read_sumstats_subgroups(subgroups, dest_dir)

#Restructure data
restructure_sumstats_subgroups(subgroups, dest_dir)

# Generate data for munge
munge_sumstats_subgroups(subgroups, dest_dir)

#Generate data including [allele frequencies] from [merged_CLEAN2.txt] ([merged_CLEAN4.txt])
metaanalysis_sumstats_subgroups(subgroups, dest_dir)


#Define subgroup
subgroups <- paste0("birthyear_bin", 9:13)

#Read summary statistics data
read_sumstats_subgroups(subgroups, dest_dir)

#Restructure data
restructure_sumstats_subgroups(subgroups, dest_dir)

# Generate data for munge
munge_sumstats_subgroups(subgroups, dest_dir)

#Generate data including [allele frequencies] from [merged_CLEAN2.txt] ([merged_CLEAN4.txt])
metaanalysis_sumstats_subgroups(subgroups, dest_dir)


#Define subgroup
subgroups <- paste0("birthyear_bin", 14:19)

#Read summary statistics data
read_sumstats_subgroups(subgroups, dest_dir)

#Restructure data
restructure_sumstats_subgroups(subgroups, dest_dir)

# Generate data for munge
munge_sumstats_subgroups(subgroups, dest_dir)

#Generate data including [allele frequencies] from [merged_CLEAN2.txt] ([merged_CLEAN4.txt])
metaanalysis_sumstats_subgroups(subgroups, dest_dir)



