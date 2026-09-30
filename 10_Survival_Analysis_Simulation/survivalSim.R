# Testing out a survival bias simulation to examine how survival bias will affect GWAS coefficients and allele frequencies.
#install.packages("meta")
library(meta)

set.seed(222)

nRep <- 1000

results <- as.data.frame(matrix(NA, nRep, 69))

for (i in 1:nRep){
N <- 1000000                                                                              # specify a sample size
maf <- .3


                                                                                          #
snp <- rbinom(N, 2, maf)                                                                   # simulate a SNP
age <- rnorm(N, 60, 10)/10                                                                # simulate age
                                                                                          #
                                                                                          # Define: 
b_0   <- .3                                                                               # Intercept
b_snp <- .025                                                                               # genetic association 
b_age <- .1                                                                               # age association
                                                                                          #
resVal <- 1- (b_snp^2*var(snp) + b_age^2*var(age))                                        # figure out the value of the residual so that the DV has a variance near 1
                                                                                          #
resid <- rnorm(N, 0, resVal)                                                              # simulate the residual
                                                                                          #
Y <- b_0 + b_snp*snp + b_age*age + resid                                                  # Simulate Y given the snp and age (no interaction)
                                                                                          #
dat1 <- data.frame(Y, snp, age)                                                           # put it into a dataframe
full <- summary(lm(Y ~ snp + age, data = dat1))                                                            # run the regression
                                                                                          #
CUTS <- quantile(age, prob = c(.2, .4, .6, .8))                                           # figure out the quintiles

no_sel <- matrix(NA, 5, 3)
# Run the regression in each quintile
                                                                                          
no_sel[1,] <- t(summary(lm(Y ~ snp + age, data = dat1[dat1$age < CUTS[1],]))$coef[,1])                         
no_sel[2,] <- t(summary(lm(Y ~ snp + age, data = dat1[dat1$age > CUTS[1] & dat1$age < CUTS[2],]))$coef[,1]) 
no_sel[3,] <- t(summary(lm(Y ~ snp + age, data = dat1[dat1$age > CUTS[2] & dat1$age < CUTS[3],]))$coef[,1]) 
no_sel[4,] <- t(summary(lm(Y ~ snp + age, data = dat1[dat1$age > CUTS[3] & dat1$age < CUTS[4],]))$coef[,1]) 
no_sel[5,] <- t(summary(lm(Y ~ snp + age, data = dat1[dat1$age > CUTS[4],]))$coef[,1])                      

no_sel

table(dat1$snp[dat1$age < CUTS[1]])
table(dat1$snp[dat1$age > CUTS[1] & dat1$age < CUTS[2]])
table(dat1$snp[dat1$age > CUTS[2] & dat1$age < CUTS[3]])
table(dat1$snp[dat1$age > CUTS[3] & dat1$age < CUTS[4]])
table(dat1$snp[dat1$age > CUTS[4]])                    

ac_noSel1 <- sum(dat1$snp[dat1$age < CUTS[1]])
ac_noSel2 <- sum(dat1$snp[dat1$age > CUTS[1] & dat1$age < CUTS[2]])
ac_noSel3 <- sum(dat1$snp[dat1$age > CUTS[2] & dat1$age < CUTS[3]])
ac_noSel4 <- sum(dat1$snp[dat1$age > CUTS[3] & dat1$age < CUTS[4]])
ac_noSel5 <- sum(dat1$snp[dat1$age > CUTS[4]])                    


n_noSel1 <- sum(table(dat1$snp[dat1$age < CUTS[1]]))
n_noSel2 <- sum(table(dat1$snp[dat1$age > CUTS[1] & dat1$age < CUTS[2]]))
n_noSel3 <- sum(table(dat1$snp[dat1$age > CUTS[2] & dat1$age < CUTS[3]]))
n_noSel4 <- sum(table(dat1$snp[dat1$age > CUTS[3] & dat1$age < CUTS[4]]))
n_noSel5 <- sum(table(dat1$snp[dat1$age > CUTS[4]])                    )

maf_noSel1 <- ac_noSel1 / (2* n_noSel1)
maf_noSel2 <- ac_noSel2 / (2* n_noSel2)
maf_noSel3 <- ac_noSel3 / (2* n_noSel3)
maf_noSel4 <- ac_noSel4 / (2* n_noSel4)
maf_noSel5 <- ac_noSel5 / (2* n_noSel5)


### Impose some sort of selection bias

b_haz  <- 1.5    # how strongly Y elevates age
lambda <- 0.05   # baseline hazard scale

H <- lambda * age * exp(b_haz * (Y - mean(Y)))    # cumulative hazard to current age
S <- exp(-H)                                      # survival probability


survived <- rbinom(N, 1, S)
dat2 <- dat1[survived == 1, ]

table(survived)                                   # note that A LOT of people didn't survive


CUTS <- quantile(age, prob = c(.2, .4, .6, .8))   # keep cuts from dat1 so strata are comparable



sel1 <- summary(lm(Y ~ snp + age, data = dat2[dat2$age < CUTS[1],]))
sel2 <- summary(lm(Y ~ snp + age, data = dat2[dat2$age > CUTS[1] & dat2$age < CUTS[2],]))
sel3 <- summary(lm(Y ~ snp + age, data = dat2[dat2$age > CUTS[2] & dat2$age < CUTS[3],]))
sel4 <- summary(lm(Y ~ snp + age, data = dat2[dat2$age > CUTS[3] & dat2$age < CUTS[4],]))
sel5 <- summary(lm(Y ~ snp + age, data = dat2[dat2$age > CUTS[4],]))


sel <- matrix(NA, 5, 3)
sel_se <- matrix(NA, 5, 3)

sel[1,] <- t(sel1$coef[,1])
sel[2,] <- t(sel2$coef[,1])
sel[3,] <- t(sel3$coef[,1])
sel[4,] <- t(sel4$coef[,1])
sel[5,] <- t(sel5$coef[,1])

sel_se[1,] <- t(sel1$coef[,2])
sel_se[2,] <- t(sel2$coef[,2])
sel_se[3,] <- t(sel3$coef[,2])
sel_se[4,] <- t(sel4$coef[,2])
sel_se[5,] <- t(sel5$coef[,2])


#sel
#sel_se

# We see a decrease (or change) in the genetic association in older age groups


mean(dat2$snp[dat2$age < CUTS[1]])
mean(dat2$snp[dat2$age > CUTS[1] & dat2$age < CUTS[2]])
mean(dat2$snp[dat2$age > CUTS[2] & dat2$age < CUTS[3]])
mean(dat2$snp[dat2$age > CUTS[3] & dat2$age < CUTS[4]])
mean(dat2$snp[dat2$age > CUTS[4]])                    

table(dat2$snp[dat2$age < CUTS[1]])
table(dat2$snp[dat2$age > CUTS[1] & dat2$age < CUTS[2]])
table(dat2$snp[dat2$age > CUTS[2] & dat2$age < CUTS[3]])
table(dat2$snp[dat2$age > CUTS[3] & dat2$age < CUTS[4]])
table(dat2$snp[dat2$age > CUTS[4]])                    

sum(table(dat2$snp[dat2$age < CUTS[1]]))
sum(table(dat2$snp[dat2$age > CUTS[1] & dat2$age < CUTS[2]]))
sum(table(dat2$snp[dat2$age > CUTS[2] & dat2$age < CUTS[3]]))
sum(table(dat2$snp[dat2$age > CUTS[3] & dat2$age < CUTS[4]]))
sum(table(dat2$snp[dat2$age > CUTS[4]])                    )


ac_Sel1 <- sum(dat2$snp[dat2$age < CUTS[1]])
ac_Sel2 <- sum(dat2$snp[dat2$age > CUTS[1] & dat2$age < CUTS[2]])
ac_Sel3 <- sum(dat2$snp[dat2$age > CUTS[2] & dat2$age < CUTS[3]])
ac_Sel4 <- sum(dat2$snp[dat2$age > CUTS[3] & dat2$age < CUTS[4]])
ac_Sel5 <- sum(dat2$snp[dat2$age > CUTS[4]])                    


n_Sel1 <- sum(table(dat2$snp[dat2$age < CUTS[1]]))
n_Sel2 <- sum(table(dat2$snp[dat2$age > CUTS[1] & dat2$age < CUTS[2]]))
n_Sel3 <- sum(table(dat2$snp[dat2$age > CUTS[2] & dat2$age < CUTS[3]]))
n_Sel4 <- sum(table(dat2$snp[dat2$age > CUTS[3] & dat2$age < CUTS[4]]))
n_Sel5 <- sum(table(dat2$snp[dat2$age > CUTS[4]])                    )

maf_Sel1 <- ac_Sel1 / (2* n_Sel1)
maf_Sel2 <- ac_Sel2 / (2* n_Sel2)
maf_Sel3 <- ac_Sel3 / (2* n_Sel3)
maf_Sel4 <- ac_Sel4 / (2* n_Sel4)
maf_Sel5 <- ac_Sel5 / (2* n_Sel5)

# We see that the allele frequency decreases (as people with the risk allele are removed from the population)




res <- metagen(TE = sel[,2],
               seTE = sel_se[,2],
               method.tau = "REML",
               method.common.ci = "IVhet")
summary(res)


results[i, ] <- c(N, maf, b_0, b_snp, b_age, t(full$coef[,1]), t(full$coef[,2]),  ## 11
t(no_sel[,1]), t(no_sel[,2]), t(no_sel[,3]),                                      ## 15
maf_noSel1, maf_noSel2, maf_noSel3,maf_noSel4,maf_noSel5,                         ## 5
t(sel[,1]), t(sel[,2]), t(sel[,3]),                                               ## 15
t(sel_se[,1]), t(sel_se[,2]), t(sel_se[,3]),                                      ## 15 
maf_Sel1,maf_Sel2,maf_Sel3,maf_Sel4,maf_Sel5,                                     ## 5 
res$Q, res$ df.Q, res$pval.Q                                                      ## 3
)

}


colnames(results) <- 
c("N", "maf", "b_0", "b_snp", "b_age", "full_b0", "full_snp", "full_age" , "full_b0_se", "full_snp_se", "full_age_se",
"no_sel_b0_1","no_sel_b0_2","no_sel_b0_3","no_sel_b0_4","no_sel_b0_5",
"no_sel_snp_1","no_sel_snp_2","no_sel_snp_3","no_sel_snp_4","no_sel_snp_5",
"no_sel_age_1","no_sel_age_2","no_sel_age_3","no_sel_age_4","no_sel_age_5",
"maf_noSel1", "maf_noSel2", "maf_noSel3","maf_noSel4","maf_noSel5",
"sel_b0_1","sel_b0_2","sel_b0_3","sel_b0_4","sel_b0_5",
"sel_snp_1","sel_snp_2","sel_snp_3","sel_snp_4","sel_snp_5",
"sel_age_1","sel_age_2","sel_age_3","sel_age_4","sel_age_5",
"sel_b0_se_1","sel_b0_se_2","sel_b0_se_3","sel_b0_se_4","sel_b0_se_5",
"sel_snp_se_1","sel_snp_se_2","sel_snp_se_3","sel_snp_se_4","sel_snp_se_5",
"sel_age_se_1","sel_age_se_2","sel_age_se_3","sel_age_se_4","sel_age_se_5",
"maf_Sel1","maf_Sel2","maf_Sel3","maf_Sel4","maf_Sel5",                                     
"Q", "df.Q", "pval.Q")  


write.table(results, "selectionSim_large.txt", row.names = F, quote = F)
write.table(results, "selectionSim_medium.txt", row.names = F, quote = F)
write.table(results, "selectionSim_small.txt", row.names = F, quote = F)


mean(results$full_b0)


require(data.table)
large <- as.data.frame(fread("selectionSim_large.txt"))
medium <- as.data.frame(fread("selectionSim_medium.txt"))
small <- as.data.frame(fread("selectionSim_small.txt"))





par(mfrow = c(3, 2))

#par(mfrow = c(1, 2))
plot(1:5, c(mean(large$sel_snp_1), mean(large$sel_snp_2), mean(large$sel_snp_3), mean(large$sel_snp_4), mean(large$sel_snp_5)), ylim = c(0.05, .12), type = "b", 
ylab = expression(beta[SNP]), col = "red", xlab = "Quintile", main = "Large Effect Size with and \nwithout Survival Bias")
points(1:5, c(mean(large$no_sel_snp_1), mean(large$no_sel_snp_2), mean(large$no_sel_snp_3), mean(large$no_sel_snp_4), mean(large$no_sel_snp_5)), type = "b", col = "blue")
legend("topright", c("No Bias", "Survival Bias"), col = c("blue", "red"), bty = "n", lwd = 1, cex = .5)

plot(density(-log10(large$pval.Q)), main = "Density of the Heterogeneity Statistic", col = "blue", lwd = 2, xlab = expression("Heterogeneity Statistic" (-log[10]("p-value"))), cex.axis = .5)
lines(c(-log10(5e-8),-log10(5e-8)), c(-.01, .1), col = "red", lwd = 3)
table(large$pval.Q < 5e-8)


#par(mfrow = c(1, 2))
plot(1:5, c(mean(medium$sel_snp_1), mean(medium$sel_snp_2), mean(medium$sel_snp_3), mean(medium$sel_snp_4), mean(medium$sel_snp_5)), ylim = c(0.01, .06), type = "b", 
ylab = expression(beta[SNP]), col = "red", xlab = "Quintile", main = "Medium Effect Size with and \nwithout Survival Bias")
points(1:5, c(mean(medium$no_sel_snp_1), mean(medium$no_sel_snp_2), mean(medium$no_sel_snp_3), mean(medium$no_sel_snp_4), mean(medium$no_sel_snp_5)), type = "b", col = "blue")

plot(density(-log10(medium$pval.Q)), main = "Density of the Heterogeneity Statistic", col = "blue", lwd = 2, 
xlab = expression("Heterogeneity Statistic" (-log[10]("p-value"))), cex.axis = .5, 
xlim = c(0, 9))
lines(c(-log10(5e-8),-log10(5e-8)), c(-.01, .1), col = "red", lwd = 3)

table(medium$pval.Q < 5e-8)

#par(mfrow = c(1, 2))
plot(1:5, c(mean(small$sel_snp_1), mean(small$sel_snp_2), mean(small$sel_snp_3), mean(small$sel_snp_4), mean(small$sel_snp_5)), ylim = c(0.01, .04), type = "b", 
ylab = expression(beta[SNP]), col = "red", xlab = "Quintile", main = "Small Effect Size with and \nwithout Survival Bias")
points(1:5, c(mean(small$no_sel_snp_1), mean(small$no_sel_snp_2), mean(small$no_sel_snp_3), mean(small$no_sel_snp_4), mean(small$no_sel_snp_5)), type = "b", col = "blue")

plot(density(-log10(small$pval.Q)), main = "Density of the Heterogeneity Statistic", col = "blue", lwd = 2, 
xlab = expression("Heterogeneity Statistic" (-log[10]("p-value"))), cex.axis = .5, 
xlim = c(0, 9))
lines(c(-log10(5e-8),-log10(5e-8)), c(-.01, .1), col = "red", lwd = 3)

table(small$pval.Q < 5e-8)







