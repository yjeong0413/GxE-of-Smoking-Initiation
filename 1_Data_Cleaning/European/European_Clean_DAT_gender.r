R.version.string

#clean environment
rm(list=ls())

# Shows loaded packages, R version, etc.
sessionInfo()   

#Define directories for "EUROPEAN"
cat("\nSet ancestry.\n")
ANCESTRY="european"

cat("\nSet directories.\n")
snp_dir=paste0("/YOUR PATH HERE/v10.2_annotation_files/",ANCESTRY,"_ac45")
sumstat_dir=paste0("/YOUR PATH HERE/",ANCESTRY)
root_dir=paste0("/YOUR PATH HERE/",ANCESTRY,"/CLEAN")
dest_dir=paste0("/YOUR PATH HERE/",ANCESTRY,"/CLEAN/GWAS_QC")
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


#Read .dat files
read_substats_subgroups <- function(subgroup, sumstat_dir) {
for (SUBGROUP in subgroup) {
cat(paste0("\n*------------------------------- ", SUBGROUP," ----------------------------------*\n"))
cat("\n------------------------------------------------------------------------------------------\n")

cat(paste0("\n1. [ever_tobacco_user_",SUBGROUP,".dat]\n"))
cat("\n")
cat("1) Read the data.\n")
sumstats <- fread(paste0(sumstat_dir,"/ever_tobacco_user_",SUBGROUP,".dat"), sep="\t", skip=1, header=F)
cat("-Assign column names.\n")
colnames(sumstats)<-c("all.data.id", "src", "log.p", "effect", "stderr", "pass", 
                      "im.num.0", "dose.b.0", "im.num.1", "dose.b.1", 
                      "AA.0", "AB.0", "BB.0", "AA.1", "AB.1", "BB.1")

    
# Create N (=the number of individuals used to compute the association statistics for each SNP)
## Genotyped case ([src]=G) & Imputed case ([src]=I)
sumstats$N_tot <-  ifelse(sumstats$src == "G", sumstats$AA.0+sumstats$AB.0+sumstats$BB.0
                                                +sumstats$AA.1+sumstats$AB.1+sumstats$BB.1,
                   ifelse(sumstats$src == "I", sumstats$im.num.0+sumstats$im.num.1, 
                          NA) #ifelse
                         ) #ifelse

cat("Return first 5 rows of [sumstats.txt]\n")
print(head(sumstats, n=5))
cat("Return last 5 rows of [sumstats.txt]\n")
print(tail(sumstats, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstats))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstats$all.data.id)))
    
# Create N_cas (=the number of cases used to compute the association statistics for each SNP)
## Genotyped case ([src]=G) & Imputed case ([src]=I)
sumstats$N_cas <-  ifelse(sumstats$src == "G", sumstats$AA.1+sumstats$AB.1+sumstats$BB.1,
                   ifelse(sumstats$src == "I", sumstats$im.num.1, 
                              NA) #ifelse
                             ) #ifelse

cat("Return first 5 rows of [sumstats.txt]\n")
print(head(sumstats, n=5))
cat("Return last 5 rows of [sumstats.txt]\n")
print(tail(sumstats, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstats))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstats$all.data.id)))
    
# Create N_con (=the number of controls used to compute the association statistics for each SNP)
## Genotyped case ([src]=G) & Imputed case ([src]=I)
sumstats$N_con <-  ifelse(sumstats$src == "G", sumstats$AA.0+sumstats$AB.0+sumstats$BB.0,
                   ifelse(sumstats$src == "I", sumstats$im.num.0, 
                              NA) #ifelse
                             ) #ifelse

cat("Return first 5 rows of [sumstats.txt]\n")
print(head(sumstats, n=5))
cat("Return last 5 rows of [sumstats.txt]\n")
print(tail(sumstats, n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(sumstats))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(sumstats$all.data.id)))
    
    
#save results
data_name <- paste0("sumstats_",SUBGROUP)
assign(data_name, sumstats, envir = .GlobalEnv)

#print the saved data
cat(paste0("Return first 5 rows of [",data_name,"]\n"))
print(head(get(data_name), n=5))
cat(paste0("Return last 5 rows of [",data_name,"]\n"))
print(tail(get(data_name), n=5))
cat("\nSize of a table: (# of rows, # of columns)\n")
print(dim(get(data_name)))
cat("\nCount if there are any missing values in 'all.data.id:\n'")
print(sum(is.na(get(data_name)$all.data.id)))

cat("\n------------------------------------------------------------------------------------------\n\n")

} #for:SUBGROUP
} #function


#Create descriptives
descriptives_subgroups <- function(subgroup) {
for (SUBGROUP in subgroup) {
#cat(paste0("\n*------------------------------- ", SUBGROUP," ----------------------------------*\n"))
cat("\n------------------------------------------------------------------------------------------\n")

data <- get(paste0("sumstats_",SUBGROUP))
 
cat("\n")    
#[src] 
## Bar plot
barplot(table(data$src),
        main = paste0(SUBGROUP, ": [src]"),
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\n",SUBGROUP,": Tabulate [src].")
print(table(data$src))

    
cat("\n")    
#[log.p] <- natural log
## histogram
hist(data$log.p, main = paste0(SUBGROUP, ": [log.p]"))
## summary
cat("Max. of [log.p]: ", max(data$log.p, na.rm = TRUE), "\n")
cat("Min. of [log.p]: ", min(data$log.p, na.rm = TRUE), "\n")
cat("Mean. of [log.p]: ", mean(data$log.p, na.rm = TRUE), "\n")
cat("Median. of [log.p]: ", median(data$log.p, na.rm = TRUE), "\n")    
    
    
cat("\n")    
#[effect] <- log OR
## histogram
hist(data$effect, main = paste0(SUBGROUP, ": [effect]"))
## summary
cat("Max. of [effect]: ", max(data$effect, na.rm = TRUE), "\n")
cat("Min. of [effect]: ", min(data$effect, na.rm = TRUE), "\n")
cat("Mean. of [effect]: ", mean(data$effect, na.rm = TRUE), "\n")
cat("Median. of [effect]: ", median(data$effect, na.rm = TRUE), "\n")


cat("\n")    
#[stderr]
## histogram
hist(data$stderr, main = paste0(SUBGROUP, ": [stderr]"))
## summary
cat("Max. of [stderr]: ", max(data$effect, na.rm = TRUE), "\n")
cat("Min. of [stderr]: ", min(data$effect, na.rm = TRUE), "\n")
cat("Mean. of [stderr]: ", mean(data$effect, na.rm = TRUE), "\n")
cat("Median. of [stderr]: ", median(data$effect, na.rm = TRUE), "\n")   
    
    
cat("\n")    
#[pass]
## Bar plot
barplot(table(data$pass), 
        main = paste0(SUBGROUP, ": [pass]"),
        col = "lightgrey", 
        xlab = "Status", 
        ylab = "Frequency")
##Tabulate
cat("\n",SUBGROUP,": Tabulate [pass].")
print(table(data$pass))

    
cat("\n")    
#[im.num.0]
## histogram
hist(data$im.num.0, main = paste0(SUBGROUP, ": [im.num.0]"))
## summary
cat("Max. of [im.num.0]: ", max(data$im.num.0, na.rm = TRUE), "\n")
cat("Min. of [im.num.0]: ", min(data$im.num.0, na.rm = TRUE), "\n")
cat("Mean. of [im.num.0]: ", mean(data$im.num.0, na.rm = TRUE), "\n")
cat("Median. of [im.num.0]: ", median(data$im.num.0, na.rm = TRUE), "\n")    

    
cat("\n")    
#[dose.b.0]
## histogram
hist(data$dose.b.0, main = paste0(SUBGROUP, ": [dose.b.0]"))
## summary
cat("Max. of [dose.b.0]: ", max(data$dose.b.0, na.rm = TRUE), "\n")
cat("Min. of [dose.b.0]: ", min(data$dose.b.0, na.rm = TRUE), "\n")
cat("Mean. of [dose.b.0]: ", mean(data$dose.b.0, na.rm = TRUE), "\n")
cat("Median. of [dose.b.0]: ", median(data$dose.b.0, na.rm = TRUE), "\n")    
    
    
cat("\n")    
#[im.num.1]
## histogram
hist(data$im.num.1, main = paste0(SUBGROUP, ": [im.num.1]"))
## summary
cat("Max. of [im.num.1]: ", max(data$im.num.1, na.rm = TRUE), "\n")
cat("Min. of [im.num.1]: ", min(data$im.num.1, na.rm = TRUE), "\n")
cat("Mean. of [im.num.1]: ", mean(data$im.num.1, na.rm = TRUE), "\n")
cat("Median. of [im.num.1]: ", median(data$im.num.1, na.rm = TRUE), "\n")
    
    
cat("\n")    
#[dose.b.1]
## histogram
#hist(data$dose.b.1, main = paste0(SUBGROUP, ": [dose.b.1]"))
## summary
cat("Max. of [dose.b.1]: ", max(data$dose.b.1, na.rm = TRUE), "\n")
cat("Min. of [dose.b.1]: ", min(data$dose.b.1, na.rm = TRUE), "\n")
cat("Mean. of [dose.b.1]: ", mean(data$dose.b.1, na.rm = TRUE), "\n")
cat("Median. of [dose.b.1]: ", median(data$dose.b.1, na.rm = TRUE), "\n")    
    
    
cat("\n")    
#[AA.0]
## histogram
hist(data$AA.0, main = paste0(SUBGROUP, ": [AA.0]"))
## summary
cat("Max. of [AA.0]: ", max(data$AA.0, na.rm = TRUE), "\n")
cat("Min. of [AA.0]: ", min(data$AA.0, na.rm = TRUE), "\n")
cat("Mean. of [AA.0]: ", mean(data$AA.0, na.rm = TRUE), "\n")
cat("Median. of [AA.0]: ", median(data$AA.0, na.rm = TRUE), "\n")
    

cat("\n")    
#[AB.0]
## histogram
hist(data$AB.0, main = paste0(SUBGROUP, ": [AB.0]"))
## summary
cat("Max. of [AB.0]: ", max(data$AB.0, na.rm = TRUE), "\n")
cat("Min. of [AB.0]: ", min(data$AB.0, na.rm = TRUE), "\n")
cat("Mean. of [AB.0]: ", mean(data$AB.0, na.rm = TRUE), "\n")
cat("Median. of [AB.0]: ", median(data$AB.0, na.rm = TRUE), "\n")
    
    
cat("\n")    
#[BB.0]
## histogram
hist(data$BB.0, main = paste0(SUBGROUP, ": [BB.0]"))
## summary
cat("Max. of [BB.0]: ", max(data$BB.0, na.rm = TRUE), "\n")
cat("Min. of [BB.0]: ", min(data$BB.0, na.rm = TRUE), "\n")
cat("Mean. of [BB.0]: ", mean(data$BB.0, na.rm = TRUE), "\n")
cat("Median. of [BB.0]: ", median(data$BB.0, na.rm = TRUE), "\n")
    
    
cat("\n")    
#[AA.1]
## histogram
hist(data$AA.1, main = paste0(SUBGROUP, ": [AA.1]"))
## summary
cat("Max. of [AA.1]: ", max(data$AA.1, na.rm = TRUE), "\n")
cat("Min. of [AA.1]: ", min(data$AA.1, na.rm = TRUE), "\n")
cat("Mean. of [AA.1]: ", mean(data$AA.1, na.rm = TRUE), "\n")
cat("Median. of [AA.1]: ", median(data$AA.1, na.rm = TRUE), "\n")
    
    
cat("\n")    
#[AB.1]
## histogram
hist(data$AB.1, main = paste0(SUBGROUP, ": [AB.1]"))
## summary
cat("Max. of [AB.1]: ", max(data$AB.1, na.rm = TRUE), "\n")
cat("Min. of [AB.1]: ", min(data$AB.1, na.rm = TRUE), "\n")
cat("Mean. of [AB.1]: ", mean(data$AB.1, na.rm = TRUE), "\n")
cat("Median. of [AB.1]: ", median(data$AB.1, na.rm = TRUE), "\n")    
    
    
cat("\n")    
#[BB.1]
## histogram
hist(data$BB.1, main = paste0(SUBGROUP, ": [BB.1]"))
## summary
cat("Max. of [BB.1]: ", max(data$BB.1, na.rm = TRUE), "\n")
cat("Min. of [BB.1]: ", min(data$BB.1, na.rm = TRUE), "\n")
cat("Mean. of [BB.1]: ", mean(data$BB.1, na.rm = TRUE), "\n")
cat("Median. of [BB.1]: ", median(data$BB.1, na.rm = TRUE), "\n")

    
cat("\n")    
#[N_tot]
temp <- data$N_tot
## histogram
hist(temp, main = paste0(SUBGROUP, ": [N_tot]"))
## summary
cat("Max. of [N_tot]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_tot]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_tot]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_tot]: ", median(temp, na.rm = TRUE), "\n")


cat("\n")    
#[N_cas]
temp <- data$N_cas
## histogram
hist(temp, main = paste0(SUBGROUP, ": [N_cas]"))
## summary
cat("Max. of [N_cas]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_cas]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_cas]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_cas]: ", median(temp, na.rm = TRUE), "\n")


cat("\n")    
#[N_con]
temp <- data$N_con
## histogram
hist(temp, main = paste0(SUBGROUP, ": [N_con]"))
## summary
cat("Max. of [N_con]: ", max(temp, na.rm = TRUE), "\n")
cat("Min. of [N_con]: ", min(temp, na.rm = TRUE), "\n")
cat("Mean. of [N_con]: ", mean(temp, na.rm = TRUE), "\n")
cat("Median. of [N_con]: ", median(temp, na.rm = TRUE), "\n")


cat("\n------------------------------------------------------------------------------------------\n\n")

} #for:SUBGROUP
} #function



#Define subsample
subgroup = c("female", "male")

#Read data
read_substats_subgroups(subgroup, sumstat_dir)

#Descriptives
descriptives_subgroups(subgroup)

# get the list of loaded data
ls()

#Extract data that has valid results
for (SUBGROUP in subgroup) {

cat(paste0("\n*------------------------------- ", SUBGROUP," ----------------------------------*\n"))
cat("\n------------------------------------------------------------------------------------------\n")
    #Define name
    data_name <- paste0("sumstat_",SUBGROUP,"_valid")
    
    #Assign data
    assign(data_name, filter(get(paste0("sumstats_",SUBGROUP)), !is.na(log.p)), envir = .GlobalEnv)
    
    #Get data
    data_valid <- get(data_name)

    cat(paste0("\nReturn first 5 rows of [",data_name,"]\n"))
    print(head(data_valid, n=5))
    cat(paste0("\nReturn last 5 rows of [",data_name,"]\n"))
    print(tail(data_valid, n=5))
    cat("\nSize of a table: (# of rows, # of columns)\n")
    print(dim(data_valid))
    cat("\nCount if there are any missing values in 'all.data.id:\n'")
    print(sum(is.na(data_valid$all.data.id)))
    
    #Save the data
    cat(paste0("Save the results in [ever_tobacco_user_",SUBGROUP,"_valid.txt].\n"))
    write.table(data_valid, file=paste0(root_dir,"/ever_tobacco_user_",SUBGROUP,"/ever_tobacco_user_",SUBGROUP,"_valid.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    cat("\n")
    cat("\n------------------------------------------------------------------------------------------\n\n")
    
}#for:


#Extract and save data passed 23andme QC using [pass]
for (SUBGROUP in subgroup) {

cat(paste0("\n*------------------------------- ", SUBGROUP," ----------------------------------*\n"))
cat("\n------------------------------------------------------------------------------------------\n")
    
    #Define name
    data_name <- paste0("sumstat_",SUBGROUP,"_passed")
    
    #Assign data
    assign(data_name, filter(get(paste0("sumstats_",SUBGROUP)), pass=="Y"), envir = .GlobalEnv)
    
    #Get data
    data_passed <- get(data_name)

    cat(paste0("\nReturn first 5 rows of [",data_name,"]\n"))
    print(head(data_passed, n=5))
    cat(paste0("\nReturn last 5 rows of [",data_name,"]\n"))
    print(tail(data_passed, n=5))
    cat("\nSize of a table: (# of rows, # of columns)\n")
    print(dim(data_passed))
    cat("\nCount if there are any missing values in 'all.data.id:\n'")
    print(sum(is.na(data_passed$all.data.id)))

    #Save the data
    cat("\nSave the results in [ever_tobacco_user_",SUBGROUP,"_passed.txt].\n")
    write.table(data_passed, file=paste0(dest_dir,"/ever_tobacco_user_",SUBGROUP,"/ever_tobacco_user_",SUBGROUP,"_passed.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    cat("\n")

cat("\n------------------------------------------------------------------------------------------\n\n")

}#for




