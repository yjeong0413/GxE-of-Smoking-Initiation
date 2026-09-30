R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#Define directories for "AFRICAN AMERICAN"
cat("Set ancestries\n")
ANCESTRY="african_american"

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


#Read dat files
cat("\n*------------------------------- Full sample files ----------------------------------*\n")
cat("\n------------------------------------------------------------------------------------------\n")

# Full Sample
cat("1. [ever_tobacco_user_passed.txt]\n")
cat("\n")
cat("1) Read the data.\n")
sumstat_full <- fread(paste0(dest_dir,"/ever_tobacco_user/ever_tobacco_user_passed.txt"), sep="\t", skip=1, header=F)
cat("-Assign column names.\n")
colnames(sumstat_full)<-c("all.data.id", "src", "log.p", "effect", "stderr", "pass", 
                          "im.num.0", "dose.b.0", "im.num.1", "dose.b.1", 
                          "AA.0", "AB.0", "BB.0", "AA.1", "AB.1", "BB.1", 
                          "N_tot", "N_cas", "N_con")

cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full$all.data.id)))

# Genotyped SNPs ([src]==G)
## Total number of A alleles (N_A)
### control
sumstat_full$N_A0 <- NA  # initialize
sumstat_full$N_A0 <- ifelse(
                            sumstat_full$src == "G",
                            ( 2*sumstat_full$AA.0+sumstat_full$AB.0),
                             sumstat_full$N_A0 # keep existing values.
                            ) #ifelse
### case
sumstat_full$N_A1 <- NA  # initialize
sumstat_full$N_A1 <- ifelse(
                            sumstat_full$src == "G",
                            ( 2*sumstat_full$AA.1+sumstat_full$AB.1),
                             sumstat_full$N_A1 # keep existing values.
                            ) #ifelse


## Total number of B alleles (N_B)
### control
sumstat_full$N_B0 <- NA  # initialize
sumstat_full$N_B0 <- ifelse(
                            sumstat_full$src == "G",
                            ( 2*sumstat_full$BB.0+sumstat_full$AB.0),
                             sumstat_full$N_B0 # keep existing values.
                            ) #ifelse
### case
sumstat_full$N_B1 <- NA  # initialize
sumstat_full$N_B1 <- ifelse(
                            sumstat_full$src == "G",
                            ( 2*sumstat_full$BB.1+sumstat_full$AB.1),
                             sumstat_full$N_B1 # keep existing values.
                            ) #ifelse


## Total number of alleles in controls (N_0=N_A0+N_B0) and cases (N_1=N_A1+N_B1)
### control
sumstat_full$N_0 <- NA  # initialize
sumstat_full$N_0 <- ifelse(
                            sumstat_full$src == "G",
                            (sumstat_full$N_A0+sumstat_full$N_B0),
                             sumstat_full$N_0 # keep existing values.
                            ) #ifelse
### case
sumstat_full$N_1 <- NA  # initialize
sumstat_full$N_1 <- ifelse(
                            sumstat_full$src == "G",
                            (sumstat_full$N_A1+sumstat_full$N_B1),
                             sumstat_full$N_1 # keep existing values.
                            ) #ifelse


## frequency of A allele ([freq.a0] and [freq.a1])
### control
sumstat_full$freq.a0 <- NA  # initialize
sumstat_full$freq.a0 <- ifelse(
                               sumstat_full$src == "G",
                              (sumstat_full$N_A0/sumstat_full$N_0),
                               sumstat_full$freq.a0 # keep existing values.
                               ) #ifelse
### case
sumstat_full$freq.a1 <- NA  # initialize
sumstat_full$freq.a1 <- ifelse(
                               sumstat_full$src == "G",
                              (sumstat_full$N_A1/sumstat_full$N_1),
                               sumstat_full$freq.a1 # keep existing values.
                               ) #ifelse


## frequency of B allele ([freq.b0] and [freq.b1])
### control
sumstat_full$freq.b0 <- NA  # initialize
sumstat_full$freq.b0 <- ifelse(
                               sumstat_full$src == "G",
                              (sumstat_full$N_B0/sumstat_full$N_0),
                               sumstat_full$freq.b0 # keep existing values.
                               ) #ifelse
### case
sumstat_full$freq.b1 <- NA  # initialize
sumstat_full$freq.b1 <- ifelse(
                               sumstat_full$src == "G",
                              (sumstat_full$N_B1/sumstat_full$N_1),
                               sumstat_full$freq.b1 # keep existing values.
                               ) #ifelse


cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'freq.a:'\n")
print(sum(is.na(sumstat_full$freq.a0)))

temp <- sumstat_full %>% filter(src == "G")
cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(temp, n=10))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(temp, n=10))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(temp))
cat("\nCount if there are any missing values in 'freq.a:'\n")
print(sum(is.na(temp$freq.a0)))

# histogram of [dose.b.0] 
cat("[dose.b.0]")
temp <- sumstat_full$dose.b.0
hist(temp)

cat("\n Max.: ", max(temp, na.rm=T))
cat("\n Min.: ", min(temp, na.rm=T),"\n")

# histogram of [dose.b.1] 
cat("[dose.b.1]")
temp <- sumstat_full$dose.b.1
hist(temp)

cat("\n Max.: ", max(temp, na.rm=T))
cat("\n Min.: ", min(temp, na.rm=T),"\n")

# Imputed SNPs ([src]==I)
## frequency of B allele (freq.b0 and freq.b1)
### control
sumstat_full$freq.b0 <- ifelse(
                               sumstat_full$src == "I",
                               sumstat_full$dose.b.0,
                               sumstat_full$freq.b0 # keep existing values.
                               ) #ifelse
### case
sumstat_full$freq.b1 <- ifelse(
                               sumstat_full$src == "I",
                               sumstat_full$dose.b.1,
                               sumstat_full$freq.b1 # keep existing values.
                               ) #ifelse


## frequency of A allele (freq.a)
### control
sumstat_full$freq.a0 <- ifelse(
                               sumstat_full$src == "I",
                              (1-sumstat_full$freq.b0),
                               sumstat_full$freq.a0 # keep existing values.
                               ) #ifelse
### case
sumstat_full$freq.a1 <- ifelse(
                               sumstat_full$src == "I",
                              (1-sumstat_full$freq.b1),
                               sumstat_full$freq.a1 # keep existing values.
                               ) #ifelse

cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'freq.a:\n'")
print(sum(is.na(sumstat_full$freq.a)))

# Genotyped SNPs ([src]==G)
## Total number of A alleles (N_A)
sumstat_full$N_A <- NA  # initialize
sumstat_full$N_A <- ifelse(
                            sumstat_full$src == "G",
                            ( 2*sumstat_full$AA.0+sumstat_full$AB.0
                             +2*sumstat_full$AA.1+sumstat_full$AB.1),
                             sumstat_full$N_A # keep existing values.
                            ) #ifelse

## Total number of B alleles (N_B)
sumstat_full$N_B <- NA  # initialize
sumstat_full$N_B <- ifelse(
                            sumstat_full$src == "G",
                            ( 2*sumstat_full$BB.0+sumstat_full$AB.0
                             +2*sumstat_full$BB.1+sumstat_full$AB.1),
                             sumstat_full$N_B # keep existing values.
                            ) #ifelse

## Total number of alleles (N_AB=N_A+N_B)
sumstat_full$N_AB <- NA  # initialize
sumstat_full$N_AB <- ifelse(
                            sumstat_full$src == "G",
                            (sumstat_full$N_A+sumstat_full$N_B),
                             sumstat_full$N_AB # keep existing values.
                            ) #ifelse

## frequency of A allele (freq.a)
sumstat_full$freq.a <- NA  # initialize
sumstat_full$freq.a <- ifelse(
                               sumstat_full$src == "G",
                              (sumstat_full$N_A/sumstat_full$N_AB),
                               sumstat_full$freq.a # keep existing values.
                               ) #ifelse

## frequency of B allele (freq.b)
sumstat_full$freq.b <- NA  # initialize
sumstat_full$freq.b <- ifelse(
                               sumstat_full$src == "G",
                              (sumstat_full$N_B/sumstat_full$N_AB),
                               sumstat_full$freq.b # keep existing values.
                               ) #ifelse

cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'freq.a:\n'")
print(sum(is.na(sumstat_full$freq.a)))

# Imputed SNPs ([src]==I)
## frequency of B allele (freq.b)
sumstat_full$freq.b <- ifelse(
                               sumstat_full$src == "I",
                              ((sumstat_full$im.num.0*sumstat_full$dose.b.0+sumstat_full$im.num.1*sumstat_full$dose.b.1)/(sumstat_full$im.num.0+sumstat_full$im.num.1)),
                               sumstat_full$freq.b # keep existing values.
                               ) #ifelse

## frequency of A allele (freq.a)
sumstat_full$freq.a <- ifelse(
                               sumstat_full$src == "I",
                              (1-sumstat_full$freq.b),
                               sumstat_full$freq.a # keep existing values.
                               ) #ifelse

cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'freq.a:\n'")
print(sum(is.na(sumstat_full$freq.a)))

# check whether [freq.a] and [freq.b] are all filled in.
cat("\nCheck whether [freq.a] and [freq.b] are all filled in.\n")
cat("-[freq.a] :",sum(is.na(sumstat_full$freq.a)), "\n")
cat("-[freq.b] :",sum(is.na(sumstat_full$freq.b)), "\n")

# [freq.a] 
data <- sumstat_full$freq.a
## histogram
hist(data)
## summary
cat("Max. of [freq.a]: ", max(data, na.rm = TRUE), "\n")
cat("Min. of [freq.a]: ", min(data, na.rm = TRUE), "\n")
cat("Mean. of [freq.a]: ", mean(data, na.rm = TRUE), "\n")
cat("Median. of [freq.a]: ", median(data, na.rm = TRUE), "\n")
cat("Length of non-empty [freq.a]: ", sum(!is.na(data)))

# [freq.b] 
data <- sumstat_full$freq.b
## histogram
hist(data)
## summary
cat("Max. of [freq.b]: ", max(data, na.rm = TRUE), "\n")
cat("Min. of [freq.b]: ", min(data, na.rm = TRUE), "\n")
cat("Mean. of [freq.b]: ", mean(data, na.rm = TRUE), "\n")
cat("Median. of [freq.b]: ", median(data, na.rm = TRUE), "\n")
cat("Length of non-empty [freq.b]: ", sum(!is.na(data)))

# [freq.a]+[freq.b] 
data.a <- sumstat_full$freq.a
data.b <- sumstat_full$freq.b
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
sumstat_full$MAF <- pmin(sumstat_full$freq.a, sumstat_full$freq.b)

cat("\nReturn first 5 rows of [ever_tobacco_user_passed.txt]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user_passed.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full$all.data.id)))

# [MAF] 
data <- sumstat_full$MAF
## histogram
hist(data)
## summary
cat("Max. of [MAF]: ", max(data, na.rm = TRUE), "\n")
cat("Min. of [MAF]: ", min(data, na.rm = TRUE), "\n")
cat("Mean. of [MAF]: ", mean(data, na.rm = TRUE), "\n")
cat("Median. of [MAF]: ", median(data, na.rm = TRUE), "\n")
cat("Length of non-empty [MAF]: ", sum(!is.na(data)))

#[src] 
## Bar plot
barplot(table(sumstat_full$src), 
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\nTabulate [src].")
table(sumstat_full$src)

#[src == G]
data_G <- sumstat_full %>%
          filter(src=="G")
#[src == I]
data_I <- sumstat_full %>%
          filter(src=="I")


#[src==G] 
## Bar plot
barplot(table(data_G$src), 
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\nTabulate [src==G].")
table(data_G$src)

#[src==I] 
## Bar plot
barplot(table(data_I$src), 
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\nTabulate [src==I].")
table(data_I$src)

#[log.p] <- natural log
## histogram
hist(sumstat_full$log.p)
## summary
cat("Max. of [log.p]: ", max(sumstat_full$log.p, na.rm = TRUE), "\n")
cat("Min. of [log.p]: ", min(sumstat_full$log.p, na.rm = TRUE), "\n")
cat("Mean. of [log.p]: ", mean(sumstat_full$log.p, na.rm = TRUE), "\n")
cat("Median. of [log.p]: ", median(sumstat_full$log.p, na.rm = TRUE), "\n")

#[effect] <- log OR
## histogram
hist(sumstat_full$effect)
## summary
cat("Max. of [effect]: ", max(sumstat_full$effect, na.rm = TRUE), "\n")
cat("Min. of [effect]: ", min(sumstat_full$effect, na.rm = TRUE), "\n")
cat("Mean. of [effect]: ", mean(sumstat_full$effect, na.rm = TRUE), "\n")
cat("Median. of [effect]: ", median(sumstat_full$effect, na.rm = TRUE), "\n")

#[stderr]
## histogram
hist(sumstat_full$stderr)
## summary
cat("Max. of [stderr]: ", max(sumstat_full$effect, na.rm = TRUE), "\n")
cat("Min. of [stderr]: ", min(sumstat_full$effect, na.rm = TRUE), "\n")
cat("Mean. of [stderr]: ", mean(sumstat_full$effect, na.rm = TRUE), "\n")
cat("Median. of [stderr]: ", median(sumstat_full$effect, na.rm = TRUE), "\n")

#[pass]
## Bar plot
barplot(table(sumstat_full$pass), 
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\nTabulate [pass].")
table(sumstat_full$pass)

#[im.num.0]
## histogram
hist(sumstat_full$im.num.0)
## summary
cat("Max. of [im.num.0]: ", max(sumstat_full$im.num.0, na.rm = TRUE), "\n")
cat("Min. of [im.num.0]: ", min(sumstat_full$im.num.0, na.rm = TRUE), "\n")
cat("Mean. of [im.num.0]: ", mean(sumstat_full$im.num.0, na.rm = TRUE), "\n")
cat("Median. of [im.num.0]: ", median(sumstat_full$im.num.0, na.rm = TRUE), "\n")

#[im.num.0] if [src=G]
## histogram
hist(data_G$im.num.0)
## summary
cat("Max. of [im.num.0]: ", max(data_G$im.num.0, na.rm = TRUE), "\n")
cat("Min. of [im.num.0]: ", min(data_G$im.num.0, na.rm = TRUE), "\n")
cat("Mean. of [im.num.0]: ", mean(data_G$im.num.0, na.rm = TRUE), "\n")
cat("Median. of [im.num.0]: ", median(data_G$im.num.0, na.rm = TRUE), "\n")

#[im.num.0] if [src=I]
## histogram
hist(data_I$im.num.0)
## summary
cat("Max. of [im.num.0]: ", max(data_I$im.num.0, na.rm = TRUE), "\n")
cat("Min. of [im.num.0]: ", min(data_I$im.num.0, na.rm = TRUE), "\n")
cat("Mean. of [im.num.0]: ", mean(data_I$im.num.0, na.rm = TRUE), "\n")
cat("Median. of [im.num.0]: ", median(data_I$im.num.0, na.rm = TRUE), "\n")

#[dose.b.0]
## histogram
hist(sumstat_full$dose.b.0)
## summary
cat("Max. of [dose.b.0]: ", max(sumstat_full$dose.b.0, na.rm = TRUE), "\n")
cat("Min. of [dose.b.0]: ", min(sumstat_full$dose.b.0, na.rm = TRUE), "\n")
cat("Mean. of [dose.b.0]: ", mean(sumstat_full$dose.b.0, na.rm = TRUE), "\n")
cat("Median. of [dose.b.0]: ", median(sumstat_full$dose.b.0, na.rm = TRUE), "\n")

#[im.num.1]
## histogram
hist(sumstat_full$im.num.1)
## summary
cat("Max. of [im.num.1]: ", max(sumstat_full$im.num.1, na.rm = TRUE), "\n")
cat("Min. of [im.num.1]: ", min(sumstat_full$im.num.1, na.rm = TRUE), "\n")
cat("Mean. of [im.num.1]: ", mean(sumstat_full$im.num.1, na.rm = TRUE), "\n")
cat("Median. of [im.num.1]: ", median(sumstat_full$im.num.1, na.rm = TRUE), "\n")

#[im.num.1] if [src=G]
## histogram
hist(data_G$im.num.1)
## summary
cat("Max. of [im.num.1]: ", max(data_G$im.num.1, na.rm = TRUE), "\n")
cat("Min. of [im.num.1]: ", min(data_G$im.num.1, na.rm = TRUE), "\n")
cat("Mean. of [im.num.1]: ", mean(data_G$im.num.1, na.rm = TRUE), "\n")
cat("Median. of [im.num.1]: ", median(data_G$im.num.1, na.rm = TRUE), "\n")

#[im.num.1] if [src=I]
## histogram
hist(data_I$im.num.1)
## summary
cat("Max. of [im.num.1]: ", max(data_I$im.num.1, na.rm = TRUE), "\n")
cat("Min. of [im.num.1]: ", min(data_I$im.num.1, na.rm = TRUE), "\n")
cat("Mean. of [im.num.1]: ", mean(data_I$im.num.1, na.rm = TRUE), "\n")
cat("Median. of [im.num.1]: ", median(data_I$im.num.1, na.rm = TRUE), "\n")

#[dose.b.1]
## histogram
hist(sumstat_full$dose.b.1)
## summary
cat("Max. of [dose.b.1]: ", max(sumstat_full$dose.b.1, na.rm = TRUE), "\n")
cat("Min. of [dose.b.1]: ", min(sumstat_full$dose.b.1, na.rm = TRUE), "\n")
cat("Mean. of [dose.b.1]: ", mean(sumstat_full$dose.b.1, na.rm = TRUE), "\n")
cat("Median. of [dose.b.1]: ", median(sumstat_full$dose.b.1, na.rm = TRUE), "\n")

#[AA.0]
## histogram
hist(sumstat_full$AA.0)
## summary
cat("Max. of [AA.0]: ", max(sumstat_full$AA.0, na.rm = TRUE), "\n")
cat("Min. of [AA.0]: ", min(sumstat_full$AA.0, na.rm = TRUE), "\n")
cat("Mean. of [AA.0]: ", mean(sumstat_full$AA.0, na.rm = TRUE), "\n")
cat("Median. of [AA.0]: ", median(sumstat_full$AA.0, na.rm = TRUE), "\n")

#[AB.0]
## histogram
hist(sumstat_full$AB.0)
## summary
cat("Max. of [AB.0]: ", max(sumstat_full$AB.0, na.rm = TRUE), "\n")
cat("Min. of [AB.0]: ", min(sumstat_full$AB.0, na.rm = TRUE), "\n")
cat("Mean. of [AB.0]: ", mean(sumstat_full$AB.0, na.rm = TRUE), "\n")
cat("Median. of [AB.0]: ", median(sumstat_full$AB.0, na.rm = TRUE), "\n")

#[BB.0]
## histogram
hist(sumstat_full$BB.0)
## summary
cat("Max. of [BB.0]: ", max(sumstat_full$BB.0, na.rm = TRUE), "\n")
cat("Min. of [BB.0]: ", min(sumstat_full$BB.0, na.rm = TRUE), "\n")
cat("Mean. of [BB.0]: ", mean(sumstat_full$BB.0, na.rm = TRUE), "\n")
cat("Median. of [BB.0]: ", median(sumstat_full$BB.0, na.rm = TRUE), "\n")

#[AA.1]
## histogram
hist(sumstat_full$AA.1)
## summary
cat("Max. of [AA.1]: ", max(sumstat_full$AA.1, na.rm = TRUE), "\n")
cat("Min. of [AA.1]: ", min(sumstat_full$AA.1, na.rm = TRUE), "\n")
cat("Mean. of [AA.1]: ", mean(sumstat_full$AA.1, na.rm = TRUE), "\n")
cat("Median. of [AA.1]: ", median(sumstat_full$AA.1, na.rm = TRUE), "\n")

#[AB.1]
## histogram
hist(sumstat_full$AB.1)
## summary
cat("Max. of [AB.1]: ", max(sumstat_full$AB.1, na.rm = TRUE), "\n")
cat("Min. of [AB.1]: ", min(sumstat_full$AB.1, na.rm = TRUE), "\n")
cat("Mean. of [AB.1]: ", mean(sumstat_full$AB.1, na.rm = TRUE), "\n")
cat("Median. of [AB.1]: ", median(sumstat_full$AB.1, na.rm = TRUE), "\n")

#[BB.1]
## histogram
hist(sumstat_full$BB.1)
## summary
cat("Max. of [BB.1]: ", max(sumstat_full$BB.1, na.rm = TRUE), "\n")
cat("Min. of [BB.1]: ", min(sumstat_full$BB.1, na.rm = TRUE), "\n")
cat("Mean. of [BB.1]: ", mean(sumstat_full$BB.1, na.rm = TRUE), "\n")
cat("Median. of [BB.1]: ", median(sumstat_full$BB.1, na.rm = TRUE), "\n")

#[N_tot]
temp <- sumstat_full$N_tot
## histogram
hist(temp)
## summary
cat("Max. of [N_tot]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_tot]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_tot]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_tot]: ", median(temp, na.rm = TRUE), "\n")

#[N_cas]
temp <- sumstat_full$N_cas
## histogram
hist(temp)
## summary
cat("Max. of [N_cas]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_cas]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_cas]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_cas]: ", median(temp, na.rm = TRUE), "\n")

#[N_con]
temp <- sumstat_full$N_con
## histogram
hist(temp)
## summary
cat("Max. of [N_con]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_con]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_con]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_con]: ", median(temp, na.rm = TRUE), "\n")

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

# keep "all.data.id" in [sumstat_full] for the passed SNPs
sumstat_full_passed_allid <- subset(sumstat_full, select=c("all.data.id"))
print(head(sumstat_full_passed_allid,n=5))
print(dim(sumstat_full_passed_allid))

#Save the data
cat("Save the results in [sumstat_full_passed_allid.txt].\n")
write.table(sumstat_full_passed_allid, file=paste0(root_dir,"/snp/Passed/sumstat_full_passed_allid.txt"), sep="\t", quote=F, row.names=F, col.names=T)
cat("\n")

# keep certain columns in [sumstat_full]
sumstat_full_subset <- subset(sumstat_full, select=c("all.data.id", "src", "log.p", "effect", "stderr", 
                                                     "N_tot", "N_cas", "N_con", "pass", 
                                                     "freq.a0", "freq.a1", "freq.b0", "freq.b1",
                                                     "freq.a", "freq.b", "MAF"))
print(head(sumstat_full_subset,n=5))
print(dim(sumstat_full_subset))

# keep certain columns in [snp_all_merged]
snp_all_merged_subset <- subset(snp_all_merged, select=c("all.data.id", "assay.name.all",
                                                         "scaffold","position","alleles", 
                                                         "im.freq.a", "im.freq.b", 
                                                         "avg.rsqr", "min.rsqr", 
                                                         "gt.rate", "gt.freq.a", "gt.freq.b"))
print(head(snp_all_merged_subset,n=5))
print(dim(snp_all_merged_subset))

# merge [sumstat_full_subset] and [snp_all_merged_subset] using "all.data.id"
merged <- merge(sumstat_full_subset, snp_all_merged_subset, by="all.data.id") # I do not specify "all=T" as I only need the matching SNP.
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
cat("\nSeparate 'Alleles' column into 'Allele 1' and 'Allele 2'.\n")
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

#Check the summary statistics of [merged_CLEAN1$gt.freq.a]+[merged_CLEAN1$gt.freq.b]
merged_CLEAN1$gt.freq.total <- ifelse(
                                      merged_CLEAN1$src == "G",
                                      (merged_CLEAN1$gt.freq.a+merged_CLEAN1$gt.freq.b),
                                      NA
                                      ) #ifelse
print(head(merged_CLEAN1, n=5))
print(tail(merged_CLEAN1, n=5))
print(dim(merged_CLEAN1))
print(sum(!is.na(merged_CLEAN1$gt.freq.total)))

# [gt.freq.total]
## histogram
hist(merged_CLEAN1$gt.freq.total)
## summary
cat("Max. of [gt.freq.total]: ", max(merged_CLEAN1$gt.freq.total, na.rm = TRUE), "\n")
cat("Min. of [gt.freq.total]: ", min(merged_CLEAN1$gt.freq.total, na.rm = TRUE), "\n")
cat("Mean. of [gt.freq.total]: ", mean(merged_CLEAN1$gt.freq.total, na.rm = TRUE), "\n")
cat("Median. of [gt.freq.total]: ", median(merged_CLEAN1$gt.freq.total, na.rm = TRUE), "\n")
cat("Length of non-empty [gt.freq.total]: ", sum(!is.na(merged_CLEAN1$gt.freq.total)))

#Initialize [MAF_old]
merged_CLEAN1$MAF_old <- NA

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

#Check the summary statistics of [merged_CLEAN1$im.freq.a]+[merged_CLEAN1$im.freq.b]
merged_CLEAN1$im.freq.total <- ifelse(
                                      merged_CLEAN1$src == "I",
                                      (merged_CLEAN1$im.freq.a+merged_CLEAN1$im.freq.b),
                                      NA
                                      ) #ifelse
print(head(merged_CLEAN1, n=5))
print(tail(merged_CLEAN1, n=5))
print(dim(merged_CLEAN1))

# [im.freq.total]
## histogram
hist(merged_CLEAN1$im.freq.total)
## summary
cat("Max. of [im.freq.total]: ", max(merged_CLEAN1$im.freq.total, na.rm = TRUE), "\n")
cat("Min. of [im.freq.total]: ", min(merged_CLEAN1$im.freq.total, na.rm = TRUE), "\n")
cat("Mean. of [im.freq.total]: ", mean(merged_CLEAN1$im.freq.total, na.rm = TRUE), "\n")
cat("Median. of [im.freq.total]: ", median(merged_CLEAN1$im.freq.total, na.rm = TRUE), "\n")
cat("Length of non-empty [im.freq.total]: ", sum(!is.na(merged_CLEAN1$im.freq.total)))

#Imputed case ([src == I])
## Normalize
#merged_CLEAN1$im.pa <- ifelse(
#                            merged_CLEAN1$src == "I",
#                            merged_CLEAN1$im.freq.a/(merged_CLEAN1$im.freq.a+merged_CLEAN1$im.freq.b),
#                            NA
#                            ) #ifelse

#print(head(merged_CLEAN1, n=5))
#print(tail(merged_CLEAN1, n=5))
#print(dim(merged_CLEAN1))



#=> HOWEVER, I checked our data. For the imputed data, freq.a+freq.b sum up to 1. Hence, we do not need to normalize. Just do "MAF = minimum of [freq.a] and [freq.b]" for the imputed data as well. Keep in mind that our case can be a special case, meaning that one should always double check whether [freq.a]+ [freq.b] sum up to 1 when they are working with imputed data.

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

# [MAF_old] <- when GENOTYPED case is considered
ind <- (merged_CLEAN1$src=="G")
## histogram
hist(merged_CLEAN1$MAF_old[ind])
## summary
cat("Max. of [MAF_old]: ", max(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Min. of [MAF_old]: ", min(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Mean. of [MAF_old]: ", mean(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Median. of [MAF_old]: ", median(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Length of non-empty [MAF_old]: ", sum(!is.na(merged_CLEAN1$MAF_old[ind])))

# [MAF_old] <- when IMPUTED case is considered
ind <- (merged_CLEAN1$src=="I")
## histogram
hist(merged_CLEAN1$MAF_old[ind])
## summary
cat("Max. of [MAF_old]: ", max(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Min. of [MAF_old]: ", min(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Mean. of [MAF_old]: ", mean(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Median. of [MAF_old]: ", median(merged_CLEAN1$MAF_old[ind], na.rm = TRUE), "\n")
cat("Length of non-empty [MAF_old]: ", sum(!is.na(merged_CLEAN1$MAF_old[ind])))

#Save the data BEFORE removing SNPs by [MAF]
#cat("Save the results in [merged_CLEAN1.txt].\n")
#write.table(merged_CLEAN1, file=paste0(dest_dir,"/ever_tobacco_user/merged_CLEAN1.txt"), sep="\t", quote=F, row.names=F, col.names=T)
#cat("\n")

#Remove SNPs with MAF < 0.01
merged_CLEAN2 <- merged_CLEAN1 %>%
                    filter(MAF>=0.01) #Automatically drops rows with MAF = NA as well (since NA >= 0.01 is NA, and filter() drops NA by default)

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
cat("Save the results in [merged_CLEAN2.txt].\n")
write.table(merged_CLEAN2, file=paste0(dest_dir,"/ever_tobacco_user/merged_CLEAN2.txt"), sep="\t", quote=F, row.names=F, col.names=T)
cat("\n")

#[src]
print(table(merged_CLEAN2$src))

#[chr]
print(table(merged_CLEAN2$chr))

#Keep 'chr1-22' only.
cat("\nKeep 'chr1-22' only.\n")
merged_CLEAN3_temp <- merged_CLEAN2 %>%
                    filter(chr %in% as.character(1:22)) # Keep chr1-22 only.

print(head(merged_CLEAN3_temp, n=5))
print(tail(merged_CLEAN3_temp, n=5))
print(dim(merged_CLEAN3_temp))

#[src]
print(table(merged_CLEAN3_temp$src))

#[chr]
print(table(merged_CLEAN3_temp$chr))

#Create raw p-value.
merged_CLEAN3_temp$pval <- exp(merged_CLEAN3_temp$log.p)

print(head(merged_CLEAN3_temp, n=5))
print(tail(merged_CLEAN3_temp, n=5))
print(dim(merged_CLEAN3_temp))

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

#Take subset of [merged_CLEAN]
merged_CLEAN3 <- merged_CLEAN3_temp %>%
                    subset(select = c("chr", "position", "snpid", "A", "B", 
                                      "log.p", "pval", "effect", "stderr", 
                                      "N_tot", "N_cas", "N_con",
                                      "MAF") # select
                          ) # subset

print(head(merged_CLEAN3, n=5))
print(tail(merged_CLEAN3, n=5))
print(dim(merged_CLEAN3))


#Save the data
cat("Save the results in [merged_CLEAN3.txt].\n")
write.table(merged_CLEAN3, file=paste0(dest_dir,"/ever_tobacco_user/merged_CLEAN3.txt"), sep="\t", quote=F, row.names=F, col.names=T)
cat("\n")

# Create [merged_CLEAN4_temp]
merged_CLEAN4_temp <- merged_CLEAN2

#[src]
print(table(merged_CLEAN4_temp$src))

#[chr]
print(table(merged_CLEAN4_temp$chr))

#Create raw p-value.
merged_CLEAN4_temp$pval <- exp(merged_CLEAN4_temp$log.p)

print(head(merged_CLEAN4_temp, n=5))
print(tail(merged_CLEAN4_temp, n=5))
print(dim(merged_CLEAN4_temp))

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

#Take subset of [merged_CLEAN]
merged_CLEAN4 <- merged_CLEAN4_temp %>%
                    subset(select = c("chr", "position", "snpid", "A", "B", 
                                      "log.p", "pval", "effect", "stderr", 
                                      "N_tot", "N_cas", "N_con",
                                      "freq.a0", "freq.a1","freq.b0", "freq.b1",
                                      "freq.a", "freq.b",
                                      "MAF") # select
                          ) # subset

print(head(merged_CLEAN4, n=5))
print(tail(merged_CLEAN4, n=5))
print(dim(merged_CLEAN4))


#Save the data
cat("Save the results in [merged_CLEAN4.txt].\n")
write.table(merged_CLEAN4, file=paste0(dest_dir,"/ever_tobacco_user/merged_CLEAN4.txt"), sep="\t", quote=F, row.names=F, col.names=T)
cat("\n")


