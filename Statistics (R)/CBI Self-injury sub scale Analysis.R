library(dplyr)
library(ggplot2)
library(car)

data <- read.csv("data_Q_CBI.csv", header=TRUE)

#checks on data
head(data)

#check & change to factors
data$Gender <- as.factor(data$Gender)
data$Group <- as.factor(data$Group)
is.factor((data$Gender))
is.factor((data$Group))

##remove 'space' after male
data$Gender <- trimws(data$Gender)
data$Gender <- factor(data$Gender, levels = c("Female", "Male"))
levels(data$Gender)

##same for group
data$Group <- trimws(data$Group)
data$Group <- factor(data$Group, levels = c("SYNGAP", "TDC"))
levels(data$Group)

###check distribution of CBI SIB_Severity 
shapiro.test(data$SIB_Severity) 


###subset data by group (SYNGAP v TDC)
data_SYNGAP <- subset(data, Group == "SYNGAP")
data_TDC <- subset(data, Group == "TDC")


#######Statistics

###SIB CB & sex
CB_SIB_sex_result <- wilcox.test(SIB_Severity ~ Gender, data = data_SYNGAP)
print(CB_SIB_sex_result)


##CBI SIB & Adaptive Behaviour 

#VineABC - sig. 
vine_result_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$VineABC, method = "spearman")
print(vine_result_SIB)


##CBI & SIB & ADHD Traits (ConIN, ConHY, CBCL (T-test & DSM)) 

spearman_SYNGAP_ConIN_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$ConIN, method = "spearman")
spearman_SYNGAP_ConHY_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$ConHY, method = "spearman")
spearman_SYNGAP_CBCL_AP_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$CBCL.AP, method = "spearman")
spearman_SYNGAP_CBCL_AD_DP_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$CBCL.AD.DP, method = "spearman")
spearman_SYNGAP_CBCL_TOT <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$CBCL.Total, method = "spearman")

# Print results
spearman_SYNGAP_ConIN_SIB
spearman_SYNGAP_ConHY_SIB
spearman_SYNGAP_CBCL_AP_SIB 
spearman_SYNGAP_CBCL_AD_DP_SIB


# Check if p-values are significant with Bonferroni correction
bonferroni_alpha <- 0.05 / 4

# Print results with adjusted significance level
cat("Spearman correlation between SIB and ConIN: p-value =", spearman_SYNGAP_ConIN_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_ConIN_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and ConHY: p-value =", spearman_SYNGAP_ConHY_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_ConHY_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and CBCL.AP: p-value =", spearman_SYNGAP_CBCL_AP_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_CBCL_AP_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and CBCL.AD.DP: p-value =", spearman_SYNGAP_CBCL_AD_DP_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_CBCL_AD_DP_SIB$p.value < bonferroni_alpha, "\n")


#####run test for CBCL & IBCL combined####

###CBI SIB & Autistic traits (SRS: Total, SCI, RRB) 
spearman_SYNGAP_SRS_T_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SrsTotal, method = "spearman")
spearman_SYNGAP_SCI_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SrsSCI, method = "spearman")
spearman_SYNGAP_RRB_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SrsRRB, method = "spearman")

# Print results
spearman_SYNGAP_SRS_T_SIB
spearman_SYNGAP_SCI_SIB
spearman_SYNGAP_RRB_SIB


# Check if p-values are significant with Bonferroni correction
bonferroni_alpha <- 0.05 / 3

# Print results with adjusted significance level
cat("Spearman correlation between SIB and SRS_T: p-value =", spearman_SYNGAP_SRS_T_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SRS_T_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and SRS SCI: p-value =", spearman_SYNGAP_SCI_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SCI_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and SRS RRB: p-value =", spearman_SYNGAP_RRB_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_RRB_SIB$p.value < bonferroni_alpha, "\n\n")



###CBI SIB & Sensory Profile scores 
spearman_SYNGAP_SEN_TOT_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SEN.TOT, method = "spearman")
spearman_SYNGAP_SEN_SK_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SEN.SK, method = "spearman")
spearman_SYNGAP_SEN_AV_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SEN.AV, method = "spearman")
spearman_SYNGAP_SEN_SEN_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SEN.SEN, method = "spearman")
spearman_SYNGAP_SEN_REG_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SEN.REG, method = "spearman")

# Print results
spearman_SYNGAP_SEN_TOT_SIB
spearman_SYNGAP_SEN_SK_SIB
spearman_SYNGAP_SEN_AV_SIB 
spearman_SYNGAP_SEN_SEN_SIB
spearman_SYNGAP_SEN_REG_SIB


# Check if p-values are significant with Bonferroni correction
bonferroni_alpha <- 0.05 / 5

# Print results with adjusted significance level
cat("Spearman correlation between SIB and SEN Total: p-value =", spearman_SYNGAP_SEN_TOT_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_TOT_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and SEN SK: p-value =", spearman_SYNGAP_SEN_SK_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_SK_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and SEN AV: p-value =", spearman_SYNGAP_SEN_AV_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_AV_SIB$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between SIB and SEN SEN: p-value =", spearman_SYNGAP_SEN_SEN_SIB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_SEN_SIB$p.value < bonferroni_alpha, "\n")

cat("Spearman correlation between SIB and SEN REG: p-value =", spearman_SYNGAP_SEN_REG$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_REG$p.value < bonferroni_alpha, "\n")



###CBI SIB & Sleep 
result_sleep_SIB <- cor.test(data$SIB_Severity, data$CSH.TSD33, method = "spearman")
print(result_sleep_SIB)



###Total CB & CBCL Int/Ext
result_CBCL_Int_SIB <- cor.test(data$SIB_Severity, data$CBCL.INT, method = "spearman")
print(result_CBCL_Int_SIB)

result_CBCL_Ext_SIB <- cor.test(data$SIB_Severity, data$CBCL.EXT, method = "spearman")
print(result_CBCL_Ext_SIB)






