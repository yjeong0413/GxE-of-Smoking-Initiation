R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#Define directories for "AFRICAN_AMERICAN"
cat("\nSet ancestry.\n")
ANCESTRY="african_american"

cat("\nSet directories.\n")
snp_dir=paste0("/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/v10.2_annotation_files/",ANCESTRY,"_ac45")
sumstat_dir=paste0("/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/",ANCESTRY)
root_dir=paste0("/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN")
dest_dir=paste0("/depot/rche/data/datasets/23andme/scratch/jeong212/smoking_GxE_interactions_study_2024/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC")
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
cat("1. [ever_tobacco_user.dat]\n")
cat("\n")
cat("1) Read the data.\n")
sumstat_full <- fread(paste0(sumstat_dir,"/ever_tobacco_user.dat"), sep="\t", skip=1, header=F)
cat("-Assign column names.\n")
colnames(sumstat_full)<-c("all.data.id", "src", "log.p", "effect", "stderr", "pass", 
                          "im.num.0", "dose.b.0", "im.num.1", "dose.b.1", 
                          "AA.0", "AB.0", "BB.0", "AA.1", "AB.1", "BB.1")
cat("\nReturn first 5 rows of [ever_tobacco_user.dat]\n")
print(head(sumstat_full, n=5))
cat("\nReturn last 5 rows of [ever_tobacco_user.dat]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full$all.data.id)))

#[src] 
## Bar plot
barplot(table(sumstat_full$src), 
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\nTabulate [src].")
table(sumstat_full$src)

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

cat("Unique [src] values: ", unique(sumstat_full$src), "\n\n")

#[im.num.0], [im.num.1], [dose.b.0], and [dose.b.1]
has_value_1 <- !is.na(sumstat_full$im.num.0)|
               !is.na(sumstat_full$im.num.1)|
               !is.na(sumstat_full$dose.b.0)|
               !is.na(sumstat_full$dose.b.1)

cat("- Count if [src]==G and ([im.num.0], [im.num.1], [dose.b.0], or [dose.b.1] is NON-missing) :", sum((sumstat_full$src=="G")&(has_value_1)), "\n")
cat("- Count if [src]==I and ([im.num.0], [im.num.1], [dose.b.0], or [dose.b.1] is NON-missing) :", sum((sumstat_full$src=="I")&(has_value_1)), "\n")

#[AA.0], [AB.0], [BB.0], [AA.1], [AB.1], and [BB.1]
has_value_2 <- !is.na(sumstat_full$AA.0)|
               !is.na(sumstat_full$AB.0)|
               !is.na(sumstat_full$BB.0)|
               !is.na(sumstat_full$AA.1)|
               !is.na(sumstat_full$AB.1)|
               !is.na(sumstat_full$BB.1)

cat("- Count if [src]==G and ([AA.0], [AB.0], [BB.0], [AA.1], [AB.1], or [BB.1] is NON-missing) :", sum((sumstat_full$src=="G")&(has_value_2)), "\n")
cat("- Count if [src]==I and ([AA.0], [AB.0], [BB.0], [AA.1], [AB.1], or [BB.1] is NON-missing) :", sum((sumstat_full$src=="I")&(has_value_2)), "\n")



#[im.num.0]
## histogram
hist(sumstat_full$im.num.0)
## summary
cat("Max. of [im.num.0]: ", max(sumstat_full$im.num.0, na.rm = TRUE), "\n")
cat("Min. of [im.num.0]: ", min(sumstat_full$im.num.0, na.rm = TRUE), "\n")
cat("Mean. of [im.num.0]: ", mean(sumstat_full$im.num.0, na.rm = TRUE), "\n")
cat("Median. of [im.num.0]: ", median(sumstat_full$im.num.0, na.rm = TRUE), "\n")

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

#[dose.b.1]
## histogram
#hist(sumstat_full$dose.b.1)
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

# Create N_tot (=the number of individuals used to compute the association statistics for each SNP)
## Genotyped case ([src]=G) & Imputed case ([src]=I)
sumstat_full$N_tot <-  ifelse(sumstat_full$src == "G", sumstat_full$AA.0+sumstat_full$AB.0+sumstat_full$BB.0
                                                        +sumstat_full$AA.1+sumstat_full$AB.1+sumstat_full$BB.1,
                       ifelse(sumstat_full$src == "I", sumstat_full$im.num.0+sumstat_full$im.num.1, 
                              NA) #ifelse
                             ) #ifelse

cat("Return first 5 rows of [sumstat_full.txt]\n")
print(head(sumstat_full, n=5))
cat("Return last 5 rows of [sumstat_full.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full$all.data.id)))

# Create N_cas (=the number of cases used to compute the association statistics for each SNP)
## Genotyped case ([src]=G) & Imputed case ([src]=I)
sumstat_full$N_cas <-  ifelse(sumstat_full$src == "G", sumstat_full$AA.1+sumstat_full$AB.1+sumstat_full$BB.1,
                       ifelse(sumstat_full$src == "I", sumstat_full$im.num.1, 
                              NA) #ifelse
                             ) #ifelse

cat("Return first 5 rows of [sumstat_full.txt]\n")
print(head(sumstat_full, n=5))
cat("Return last 5 rows of [sumstat_full.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full$all.data.id)))

# Create N_con (=the number of controls used to compute the association statistics for each SNP)
## Genotyped case ([src]=G) & Imputed case ([src]=I)
sumstat_full$N_con <-  ifelse(sumstat_full$src == "G", sumstat_full$AA.0+sumstat_full$AB.0+sumstat_full$BB.0,
                       ifelse(sumstat_full$src == "I", sumstat_full$im.num.0, 
                              NA) #ifelse
                             ) #ifelse

cat("Return first 5 rows of [sumstat_full.txt]\n")
print(head(sumstat_full, n=5))
cat("Return last 5 rows of [sumstat_full.txt]\n")
print(tail(sumstat_full, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full$all.data.id)))

#[N_tot]
temp <- sumstat_full$N_tot
## histogram
hist(temp, main = paste0("full : [N_tot]"))
## summary
cat("Max. of [N_tot]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_tot]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_tot]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_tot]: ", median(temp, na.rm = TRUE), "\n")

#[N_cas]
temp <- sumstat_full$N_cas
## histogram
hist(temp, main = paste0("full : [N_cas]"))
## summary
cat("Max. of [N_cas]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_cas]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_cas]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_cas]: ", median(temp, na.rm = TRUE), "\n")

#[N_con]
temp <- sumstat_full$N_con
## histogram
hist(temp, main = paste0("full : [N_con]"))
## summary
cat("Max. of [N_con]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_con]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_con]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_con]: ", median(temp, na.rm = TRUE), "\n")

#Extract data that has valid results
sumstat_full_valid <- filter(sumstat_full, !is.na(log.p))

cat("Return first 5 rows of [sumstat_full_valid.txt]\n")
print(head(sumstat_full_valid, n=5))
cat("Return last 5 rows of [sumstat_full_valid.txt]\n")
print(tail(sumstat_full_valid, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full_valid))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full_valid$all.data.id)))

#Save the data
cat("Save the results in [ever_tobacco_user_valid.txt].\n")
write.table(sumstat_full_valid, file=paste0(root_dir,"/ever_tobacco_user/ever_tobacco_user_valid.txt"), sep="\t", quote=F, row.names=F, col.names=T)
cat("\n")



#Extract and save data passed 23andme QC using [pass]
sumstat_full_passed <- filter(sumstat_full, pass=="Y")

print(head(sumstat_full_passed, n=5))
cat("Return last 5 rows of [ever_tobacco_user.dat]\n")
print(tail(sumstat_full_passed, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstat_full_passed))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstat_full_passed$all.data.id)))

#Save the data
cat("Save the results in [ever_tobacco_user_passed.txt].\n")
write.table(sumstat_full_passed, file=paste0(dest_dir,"/ever_tobacco_user/ever_tobacco_user_passed.txt"), sep="\t", quote=F, row.names=F, col.names=T)
cat("\n")








