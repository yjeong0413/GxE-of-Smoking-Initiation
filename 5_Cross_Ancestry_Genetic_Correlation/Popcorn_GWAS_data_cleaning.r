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
                 AFR = "african_american",
                 AMR = "latino"
                )

print(names(ancestry_map))
print(ancestry_map[["EUR"]])
print(ancestry_map[["EAS"]])
print(ancestry_map[["AFR"]])
print(ancestry_map[["AMR"]])

for (ANC in names(ancestry_map)) {
  
  ANCESTRY <- ancestry_map[[ANC]]
    
  cat("\n")
  cat(paste0("\n*===================================== ", ANCESTRY, " ================================*\n"))

    # Read reference panel data
    file_path <- paste0(reference_dir,"/",ANC,"/",ANC,".chr1_22.bim")
    reference <- fread(file_path, sep = "\t", header = FALSE)

    colnames(reference) <- c("chr","ref_snpid","cm","pos","ref_a1","ref_a2")

    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(reference, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(reference, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(reference))
    
    # clean
    reference_data <- reference %>%
                        transmute(chr = as.character(chr),
                                  pos = as.integer(pos),
                                  ref_snpid,
                                  ref_a1 = toupper(ref_a1),
                                  ref_a2 = toupper(ref_a2)
                                  )

    # Reassign
    data_name <- paste0("reference_",ANC)
    assign(data_name, reference_data, envir=.GlobalEnv)
      
    cat(paste0("\nReturn first 5 rows of [", data_name, "]\n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows of [", data_name, "]\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))
    
   cat("\n==========================================================================================\n")
   cat("\n")

} #for: ANC

# Define subgroup
type <- "gender"
subgroups <- c("female", "male")
subgroups

for (ANC in names(ancestry_map)) {

  ANCESTRY <- ancestry_map[[ANC]]

  cat("\n")
  cat(paste0("\n*===================================== ", ANCESTRY, " ================================*\n"))
  for (SUBGROUP in subgroups) {
    cat("\n")
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))
    
    # Read data
    file_path <- paste0(start_dir,"/WORKPLACE/ldsc/",
                        ANCESTRY,"/CLEAN/GWAS_QC/ever_tobacco_user_", 
                        SUBGROUP,"/merged_CLEAN4.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("1) Read the data.\n")
    
    merged_data <- fread(file_path, sep = "\t", header = T)
    
    cat("\n[merged_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(merged_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(merged_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(merged_data))

    # Read hg19 positions
    file_path=paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Sumstats/",
                     ANCESTRY,"/",type,"/",ANCESTRY,"_",SUBGROUP,"_hg19.txt")
    cat(paste0("\n2. [", file_path, "]\n"))
    cat("2) Read hg19 position data.\n")
    
    hg19_data <- fread(file_path, sep = "\t", header = F)
    colnames(hg19_data) <- c("chr_hg19", "pos_start_hg19", "pos_end_hg19", "rsid_hg19")
    
    cat("\n[hg19_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(hg19_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(hg19_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(hg19_data))
      
    # Clean data
    cleaned_data <- merged_data %>%
                      inner_join(hg19_data, by=c("snpid"="rsid_hg19"))%>%
                      select(chr, pos_end_hg19, snpid, 
                             B, A, freq.a, 
                             N_tot, effect, stderr) %>%
                      rename(pos=pos_end_hg19,
                             rsid=snpid, 
                             a1=B, a2=A, af=freq.a, 
                             N=N_tot,
                             beta=effect, SE=stderr)

    cat("\n[cleaned_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(cleaned_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(cleaned_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(cleaned_data))

    # Harmonize data
     ## merge with reference panel by chr and position 
     harmonized_temp <- inner_join(cleaned_data, 
                                  get(paste0("reference_",ANC)),
                                  by = c("chr","pos"))
     ## check allele orientation
      comp <- c("A"="T", "T"="A", "C"="G", "G"="C")
      harmonized_temp <- harmonized_temp %>%
                          mutate(a1c = comp[a1],
                                 a2c = comp[a2],
                                 match_type = case_when(
                                     a1 == ref_a1 & a2 == ref_a2 ~ "match",
                                     a1 == ref_a2 & a2 == ref_a1 ~ "reverse",
                                     a1c == ref_a1 & a2c == ref_a2 ~ "strand_match",
                                     a1c == ref_a2 & a2c == ref_a1 ~ "strand_reverse",
                                     TRUE ~ "mismatch"
                                 ) #case_when
                                )#mutate
      ## harmonize beta and af
       harmonized_data <- harmonized_temp %>%
                              mutate(beta_hm = case_when(
                                      match_type %in% c("match", "strand_match") ~ beta,
                                      match_type %in% c("reverse", "strand_reverse") ~ -beta,
                                      TRUE ~ NA_real_
                                    ),
                                    # Popcorn AF should be frequency of A2.
                                    # After forcing final a2 = ref_a2:
                                    af_hm = case_when(
                                      match_type %in% c("match", "strand_match") ~ af,
                                      match_type %in% c("reverse", "strand_reverse") ~ 1 - af,
                                      TRUE ~ NA_real_
                                    ),
                                    # alleles alinged to reference panel
                                    a1_hm = ref_a1,
                                    a2_hm = ref_a2
                                    ) %>%
                              filter(!is.na(beta_hm),
                                     !is.na(af_hm),
                                     af_hm > 0,
                                     af_hm < 1)
      
      # final cleaning
      cleaned_harmonized_data <- harmonized_data %>%
                                  transmute(chr, pos,
                                            rsid=ref_snpid,
                                            a1=ref_a1,
                                            a2=ref_a2,
                                            af=af_hm,
                                            N,
                                            beta = beta_hm,
                                            SE)
    # assign data
    data_name <- paste0("cleaned_harmonized_", SUBGROUP)
    assign(data_name, cleaned_harmonized_data, envir=.GlobalEnv)
    
    cat("\n[", data_name, "]: \n")
    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

    #Save the data
    cat("Save the results in [",data_name,"].\n")
    file_path=paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Sumstats/",ANCESTRY,"/",type,"/",ANCESTRY,"_",SUBGROUP,".txt")
    write.table(get(data_name), 
                file=file_path,
                sep="\t", quote=F, row.names=F, col.names=T)

   cat("\n------------------------------------------------------------------------------------------\n")
   cat("\n")
  } #for: SUBGROUP
   cat("\n==========================================================================================\n")
   cat("\n")
} #for: ANC

# Define subgroup
type <- "birthyear"
subgroups <- paste0("birthyear_bin",4:19)
subgroups

for (ANC in names(ancestry_map)) {

  ANCESTRY <- ancestry_map[[ANC]]

  cat("\n")
  cat(paste0("\n*===================================== ", ANCESTRY, " ================================*\n"))
  for (SUBGROUP in subgroups) {
    cat("\n")
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))
    
    # Read data
    file_path <- paste0(start_dir,"/WORKPLACE/ldsc/",
                        ANCESTRY,"/CLEAN/GWAS_QC/ever_tobacco_user_", 
                        SUBGROUP,"/merged_CLEAN4.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("1) Read the data.\n")
    
    merged_data <- fread(file_path, sep = "\t", header = T)
    
    cat("\n[merged_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(merged_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(merged_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(merged_data))

    # Read hg19 positions
    file_path=paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Sumstats/",
                     ANCESTRY,"/",type,"/",ANCESTRY,"_",SUBGROUP,"_hg19.txt")
    cat(paste0("\n2. [", file_path, "]\n"))
    cat("2) Read hg19 position data.\n")
    
    hg19_data <- fread(file_path, sep = "\t", header = F)
    colnames(hg19_data) <- c("chr_hg19", "pos_start_hg19", "pos_end_hg19", "rsid_hg19")
    
    cat("\n[hg19_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(hg19_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(hg19_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(hg19_data))
      
    # Clean data
    cleaned_data <- merged_data %>%
                      inner_join(hg19_data, by=c("snpid"="rsid_hg19"))%>%
                      select(chr, pos_end_hg19, snpid, 
                             B, A, freq.a, 
                             N_tot, effect, stderr) %>%
                      rename(pos=pos_end_hg19,
                             rsid=snpid, 
                             a1=B, a2=A, af=freq.a, 
                             N=N_tot,
                             beta=effect, SE=stderr)

    cat("\n[cleaned_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(cleaned_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(cleaned_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(cleaned_data))

    # Harmonize data
     ## merge with reference panel by chr and position 
     harmonized_temp <- inner_join(cleaned_data, 
                                  get(paste0("reference_",ANC)),
                                  by = c("chr","pos"))
     ## check allele orientation
      comp <- c("A"="T", "T"="A", "C"="G", "G"="C")
      harmonized_temp <- harmonized_temp %>%
                          mutate(a1c = comp[a1],
                                 a2c = comp[a2],
                                 match_type = case_when(
                                     a1 == ref_a1 & a2 == ref_a2 ~ "match",
                                     a1 == ref_a2 & a2 == ref_a1 ~ "reverse",
                                     a1c == ref_a1 & a2c == ref_a2 ~ "strand_match",
                                     a1c == ref_a2 & a2c == ref_a1 ~ "strand_reverse",
                                     TRUE ~ "mismatch"
                                 ) #case_when
                                )#mutate
      ## harmonize beta and af
       harmonized_data <- harmonized_temp %>%
                              mutate(beta_hm = case_when(
                                      match_type %in% c("match", "strand_match") ~ beta,
                                      match_type %in% c("reverse", "strand_reverse") ~ -beta,
                                      TRUE ~ NA_real_
                                    ),
                                    # Popcorn AF should be frequency of A2.
                                    # After forcing final a2 = ref_a2:
                                    af_hm = case_when(
                                      match_type %in% c("match", "strand_match") ~ af,
                                      match_type %in% c("reverse", "strand_reverse") ~ 1 - af,
                                      TRUE ~ NA_real_
                                    ),
                                    # alleles alinged to reference panel
                                    a1_hm = ref_a1,
                                    a2_hm = ref_a2
                                    ) %>%
                              filter(!is.na(beta_hm),
                                     !is.na(af_hm),
                                     af_hm > 0,
                                     af_hm < 1)
      
      # final cleaning
      cleaned_harmonized_data <- harmonized_data %>%
                                  transmute(chr, pos,
                                            rsid=ref_snpid,
                                            a1=ref_a1,
                                            a2=ref_a2,
                                            af=af_hm,
                                            N,
                                            beta = beta_hm,
                                            SE)
    # assign data
    data_name <- paste0("cleaned_harmonized_", SUBGROUP)
    assign(data_name, cleaned_harmonized_data, envir=.GlobalEnv)
    
    cat("\n[", data_name, "]: \n")
    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

    #Save the data
    cat("Save the results in [",data_name,"].\n")
    file_path=paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Sumstats/",ANCESTRY,"/",type,"/",ANCESTRY,"_",SUBGROUP,".txt")
    write.table(get(data_name), 
                file=file_path,
                sep="\t", quote=F, row.names=F, col.names=T)

   cat("\n------------------------------------------------------------------------------------------\n")
   cat("\n")
  } #for: SUBGROUP
   cat("\n==========================================================================================\n")
   cat("\n")
} #for: ANC

# Define subgroup
type <- "birthyear"
subgroups <- paste0("groupingBIN",1:4)
subgroups

for (ANC in names(ancestry_map)) {

  ANCESTRY <- ancestry_map[[ANC]]

  cat("\n")
  cat(paste0("\n*===================================== ", ANCESTRY, " ================================*\n"))
  for (SUBGROUP in subgroups) {
    cat("\n")
    cat(paste0("\n*----------------------------------- ", SUBGROUP, " --------------------------------*\n"))

    # Read data
    MA_dir=paste0(start_dir,"/WORKPLACE/ldsc/",ANCESTRY,"/CLEAN/GWAS_QC/Meta_Analysis")
    file_path <- paste0(MA_dir,"/",type,"/METAL_",ANCESTRY,"_",type,"_MAresults_",SUBGROUP,"_cleaned.txt")
    cat(paste0("\n1. [", file_path, "]\n"))
    cat("1) Read the data.\n")
    
    merged_data <- fread(file_path, sep = "\t", header = T)
    
    cat("\n[merged_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(merged_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(merged_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(merged_data))

    # Read hg19 positions
    file_path=paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Sumstats/",
                     ANCESTRY,"/",type,"/",ANCESTRY,"_",SUBGROUP,"_hg19.txt")
    cat(paste0("\n2. [", file_path, "]\n"))
    cat("2) Read hg19 position data.\n")
    
    hg19_data <- fread(file_path, sep = "\t", header = F)
    colnames(hg19_data) <- c("chr_hg19", "pos_start_hg19", "pos_end_hg19", "rsid_hg19")
    
    cat("\n[hg19_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(hg19_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(hg19_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(hg19_data))
      
    # Clean data
    cleaned_data <- merged_data %>%
                      inner_join(hg19_data, by=c("MarkerName"="rsid_hg19"))%>%
                      select(chr, pos_end_hg19, MarkerName, 
                             B, A, Freq1_new, 
                             N_tot, Effect_new, StdErr) %>%
                      mutate(freq.a = 1-Freq1_new) %>%
                      rename(pos=pos_end_hg19,
                             rsid=MarkerName, 
                             a1=B, a2=A, af=freq.a, 
                             N=N_tot,
                             beta=Effect_new, SE=StdErr) %>%
                      transmute(chr = as.character(chr),
                                pos = as.integer(pos),
                                rsid,
                                a1 = toupper(a1),
                                a2 = toupper(a2),
                                af, N, beta, SE)

    cat("\n[cleaned_data]: \n")
    cat(paste0("\nReturn first 5 rows:\n"))
    print(head(cleaned_data, n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(cleaned_data, n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(cleaned_data))

    # Harmonize data
     ## merge with reference panel by chr and position 
     harmonized_temp <- inner_join(cleaned_data, 
                                  get(paste0("reference_",ANC)),
                                  by = c("chr","pos"))
     ## check allele orientation
      comp <- c("A"="T", "T"="A", "C"="G", "G"="C")
      harmonized_temp <- harmonized_temp %>%
                          mutate(a1c = comp[a1],
                                 a2c = comp[a2],
                                 match_type = case_when(
                                     a1 == ref_a1 & a2 == ref_a2 ~ "match",
                                     a1 == ref_a2 & a2 == ref_a1 ~ "reverse",
                                     a1c == ref_a1 & a2c == ref_a2 ~ "strand_match",
                                     a1c == ref_a2 & a2c == ref_a1 ~ "strand_reverse",
                                     TRUE ~ "mismatch"
                                 ) #case_when
                                )#mutate
      ## harmonize beta and af
       harmonized_data <- harmonized_temp %>%
                              mutate(beta_hm = case_when(
                                      match_type %in% c("match", "strand_match") ~ beta,
                                      match_type %in% c("reverse", "strand_reverse") ~ -beta,
                                      TRUE ~ NA_real_
                                    ),
                                    # Popcorn AF should be frequency of A2.
                                    # After forcing final a2 = ref_a2:
                                    af_hm = case_when(
                                      match_type %in% c("match", "strand_match") ~ af,
                                      match_type %in% c("reverse", "strand_reverse") ~ 1 - af,
                                      TRUE ~ NA_real_
                                    ),
                                    # alleles alinged to reference panel
                                    a1_hm = ref_a1,
                                    a2_hm = ref_a2
                                    ) %>%
                              filter(!is.na(beta_hm),
                                     !is.na(af_hm),
                                     af_hm > 0,
                                     af_hm < 1)
      
      # final cleaning
      cleaned_harmonized_data <- harmonized_data %>%
                                  transmute(chr, pos,
                                            rsid=ref_snpid,
                                            a1=ref_a1,
                                            a2=ref_a2,
                                            af=af_hm,
                                            N,
                                            beta = beta_hm,
                                            SE)
    # assign data
    data_name <- paste0("cleaned_harmonized_", SUBGROUP)
    assign(data_name, cleaned_harmonized_data, envir=.GlobalEnv)
    
    cat("\n[", data_name, "]: \n")
    cat(paste0("\nReturn first 5 rows: \n"))
    print(head(get(data_name), n = 5))    
    cat(paste0("\nReturn last 5 rows:\n"))
    print(tail(get(data_name), n = 5))
    cat("\nSize of table: (# of rows, # of columns)\n")
    print(dim(get(data_name)))

    #Save the data
    cat("Save the results in [",data_name,"].\n")
    file_path=paste0(start_dir,"/WORKPLACE/ldsc/Popcorn/Sumstats/",ANCESTRY,"/",type,"/",ANCESTRY,"_",SUBGROUP,".txt")
    write.table(get(data_name), 
                file=file_path,
                sep="\t", quote=F, row.names=F, col.names=T)
      
   cat("\n------------------------------------------------------------------------------------------\n")
   cat("\n")
  } #for: SUBGROUP
   cat("\n==========================================================================================\n")
   cat("\n")
} #for: ANC


