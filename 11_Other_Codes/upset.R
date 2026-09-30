## Upset Plot for SI GxE results

# load the required packages
require(data.table)
#require(fuzzyjoin)
#require(dplyr)
require(UpSetR)


res <- as.data.frame(fread("locus_level_highPIP_patterns_core_priority.tsv"))
#alt <- as.data.frame(fread("locus_level_highPIP_patterns_core_membership.tsv"))

#We didn't mean to finemap the generation / sex loci where there wasn't a significant signal

table(res$is_generation_moderated, useNA = "ifany")
table(res$is_sex_moderated, useNA = "ifany")
	
table(res$is_generation_moderated	 ==F & res$is_generation_main ==F, res$generation_highPIP_pattern_core)
table(res$is_sex_moderated ==F & res$is_sex_main ==F, res$sex_highPIP_pattern_core)

res$generation_highPIP_pattern_core[res$is_generation_moderated ==F & res$is_generation_main ==F] <- NA
res$generation_highPIP_pattern_core[is.na(res$is_generation_moderated) & is.na(res$is_generation_main)] <- NA

res$sex_highPIP_pattern_core[res$is_sex_moderated ==F & res$is_sex_main ==F] <- NA
res$sex_highPIP_pattern_core[is.na(res$is_sex_moderated) & is.na(res$is_sex_moderated)] <- NA



table(res$is_generation_moderated ==F & res$is_generation_main ==F, res$generation_highPIP_pattern_core)
table(res$is_sex_moderated ==F & res$is_sex_main ==F, res$sex_highPIP_pattern_core)

table(alt$generation_highPIP_pattern_core)
table(alt$sex_highPIP_pattern_core)


sex_het       <- res$meta_locus_id[res$is_sex_moderated]
sex_main      <- res$meta_locus_id[res$is_sex_main]
sex_sharedPIP <- res$meta_locus_id[res$sex_highPIP_pattern_core =="Shared high-PIP"]
sex_specPIP   <- res$meta_locus_id[res$sex_highPIP_pattern_core =="Female high-PIP only"]
sex_specPIP   <- res$meta_locus_id[res$sex_highPIP_pattern_core =="Male high-PIP only"]
sex_noPIP     <- res$meta_locus_id[res$sex_highPIP_pattern_core =="No high-PIP"]


gen_het        <- res$meta_locus_id[res$is_generation_moderated]
gen_main       <- res$meta_locus_id[res$is_generation_main]
gen_sharedPIP  <- res$meta_locus_id[res$generation_highPIP_pattern_core=="Shared across generations"]
gen_specPIP    <- res$meta_locus_id[res$generation_highPIP_pattern_core=="G1_Silent only"]
gen_specPIP    <- res$meta_locus_id[res$generation_highPIP_pattern_core=="G2_BabyBoomer only"]
gen_specPIP    <- res$meta_locus_id[res$generation_highPIP_pattern_core=="G3_GenX only"]
gen_specPIP    <- res$meta_locus_id[res$generation_highPIP_pattern_core=="G4_Millennial only"]
gen_noPIP      <- res$meta_locus_id[res$generation_highPIP_pattern_core=="No high-PIP"]


length(sex_main) + length(sex_het)
length(gen_main) + length(gen_het)

sex_sharedPIP <- sex_sharedPIP [!is.na(sex_sharedPIP )]
sex_specPIP   <- sex_specPIP   [!is.na(sex_specPIP   )]

gen_sharedPIP <- gen_sharedPIP [!is.na(gen_sharedPIP )]
gen_specPIP   <- gen_specPIP   [!is.na(gen_specPIP   )]


mySNPSets <- list(
  "Sex Moderated"         = sex_het,
  "Sex Main Effect"       = sex_main,
  "Sex Shared High PIP"   = sex_sharedPIP,
  "Sex Specific High PIP" = sex_specPIP, 
  "Sex No Hight PIP"      = sex_noPIP, 

  "Generation Moderated"           = gen_het,
  "Generation Main Effect"          = gen_main,
  "Generation Shared High PIP"     = gen_sharedPIP,
  "Generation Specific High PIP"       = gen_specPIP, 
  "Generation No Hight PIP"         = gen_noPIP
  
)

# Upset plot - Just sex
upset(fromList(mySNPSets), 
      sets = rev(c("Sex Moderated", "Sex Main Effect", "Sex Shared High PIP", "Sex Specific High PIP")), 
      keep.order = TRUE, 
	  main.bar.color = c("blue", "red", "navy", "pink", "red", "maroon"), 
	  sets.bar.color = c("red", "blue", "blue", "red"), 
	  mainbar.y.label = "Intersection of Main, Moderated \nand Fine Mapping Effects \non Smoking Initiation for Sex", 
	  text.scale = 1.5)

# Upset plot - Just generation  
upset(fromList(mySNPSets), 
      sets = rev(c("Generation Moderated", "Generation Main Effect", "Generation Shared High PIP", "Generation Specific High PIP")), 
      keep.order = TRUE, 
	  main.bar.color = c("blue", "red", "navy", "pink", "red"),
	  sets.bar.color = c("red", "blue", "blue", "red"),
	  mainbar.y.label = "Intersection of Main, Moderated \nand Fine Mapping Effects \non Smoking Initiation for Generation", 
	  text.scale = 1.5)


# Upset plot - Just hetP  
upset(fromList(mySNPSets), 
      sets = rev(c("Sex Moderated", "Sex Main Effect", "Generation Moderated", "Generation Main Effect")), 
      keep.order = TRUE, 
	  main.bar.color = c("blue", "blue", "red", "navy","pink","pink", "maroon"),
	  sets.bar.color = c("red", "blue", "red", "blue"),
	  mainbar.y.label = "Intersection of Main and Moderated \nEffects on Smoking Initiation \nfor Heterogeneity across Moderators", 
	  text.scale = 1.5)


# Upset plot - Just finemap  
upset(fromList(mySNPSets), 
      sets = rev(c("Sex Shared High PIP", "Sex Specific High PIP", "Generation Shared High PIP", "Generation Specific High PIP")), 
      keep.order = TRUE, 
	  main.bar.color = c("red", "blue", "blue", "red", "navy", "pink", "maroon", "pink"),
	  sets.bar.color = c("red", "blue", "red", "blue"),
	  mainbar.y.label = "Intersection of Fine Mapping \nEffects on Smoking Initiation \nacross Moderators", 
	  text.scale = 1.5)






col2 <- 
c("blue", "blue", "red", "navy",
"red",    "pink", "navy","navy","navy",
"pink","red","red","red","red",
"red","green","red","pink","pink",
"pink","red","red","red","red",
"red")



upset(fromList(mySNPSets), 
      sets = rev(c("Sex Main Effect", "Sex Moderated", "Sex Shared High PIP",  "Sex Specific High PIP",
			 "Generation Main Effect", "Generation Moderated",  "Generation Shared High PIP","Generation Specific High PIP")), 
      keep.order = TRUE, 
	  main.bar.color = col2,
	  sets.bar.color = c("blue", "red", "blue", "red", "blue", "red", "blue", "red"), 
	  mainbar.y.label = "Intersection of Main, Moderated and \nFine Mapping Effects on Smoking Initiation \nfor Sex and Generation", 
	  text.scale = 1.5)








