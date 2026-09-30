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

anc_pairs <- combn(unname(ancestry_map), 2, simplify = FALSE)
anc_pairs

for (pair in anc_pairs){

    anc1 <- pair[1]
    anc2 <- pair[2]
    
    ANC <- paste0(anc1,"_",anc2)
    
  cat("\n")
  cat(paste0("\n*===================================== ", ANC, " ================================*\n"))
    
    # Read score data
    file_path <- paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Output/Step1_Scores/",
                        ref,"/",ANC,".cscore")
    score <- fread(file_path, sep = "\t", header = FALSE)

    colnames(score) <- c("chr","pos","rsid",
                             "a1","a2","a1f_pop1", "a1f_pop2",
                             "ld_pop1","ld_pop2",
                             "cc_score")
    cat("\n[score]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(score, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(score, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(score))

    # clean data
    score_data <- score %>%
                    transmute(
                        chr = as.character(as.integer(chr)),
                        pos = as.integer(pos),
                        rsid = paste(as.integer(chr),as.integer(pos),a2,a1, sep=":"),
                        a1 = as.character(a1),
                        a2 = as.character(a2), 
                        a1f_pop1 = as.numeric(a1f_pop1), 
                        a1f_pop2 = as.numeric(a1f_pop2),
                        ld_pop1 = as.numeric(ld_pop1), 
                        ld_pop2 = as.numeric(ld_pop2),
                        cc_score = as.numeric(cc_score)
                    )

    cat("\n[score_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(score_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(score_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(score_data))

    #Save the data
    data_name <- paste0("score_data")
    cat("Save the results in [",data_name,"].\n")
    file_path <- paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Output/Step1_Scores/",
                        ref,"/",ANC,"_cleaned.cscore")
    write.table(get(data_name), 
                file=file_path,
                sep="\t", quote=F, row.names=F, col.names=F)
    
   cat("\n==========================================================================================\n")
   cat("\n")
} #for: pair


