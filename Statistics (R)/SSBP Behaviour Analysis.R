library(dplyr)
library(ggplot2)
library(car)
library(reshape2)
library(corrplot)
library(GGally)
library(ggpubr)
library(gridExtra)
library(writexl)

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

###check distribution & do bonferroni correction for multiple comparisons

# Run Shapiro-Wilk tests and collect p-values
#NOTE - won't run bc mean_severity isn't calculated until later in script
p_values <- c(
  shapiro.test(data$Age.Months)$p.value, 
  shapiro.test(data$Total.CBI.Score)$p.value, #non-parametric
  shapiro.test(data$VineABC)$p.value, 
  shapiro.test(data$ConIN)$p.value, #non-parametric
  shapiro.test(data$ConHY)$p.value,
  shapiro.test(data$CBCL.AP)$p.value, #t-score
  shapiro.test(data$CBCL.AD.DP)$p.value, #DSM
  shapiro.test(data$IBCL.ATT)$p.value, #t-score
  shapiro.test(data$IBCL.AD.DP)$p.value, #DSM
  shapiro.test(data$CSH.TSD33)$p.value, 
  shapiro.test(data$SrsTotal)$p.value, #non-parametric
  shapiro.test(data$SrsSCI)$p.value, #non-parametric
  shapiro.test(data$SrsRRB)$p.value, #non-parametric
  shapiro.test(data$SEN.TOT)$p.value, #total sensory score 
  shapiro.test(data$SEN.SK)$p.value,
  shapiro.test(data$SEN.AV)$p.value,
  shapiro.test(data$SEN.SEN)$p.value,
  shapiro.test(data$SEN.REG)$p.value, #non-parametric 
  shapiro.test(data$CBCL.INT)$p.value,
  shapiro.test(data$CBCL.EXT)$p.value,
  shapiro.test(data$mean_severity)$p.value
)

#number of tests
n_tests <- length(p_values)

# Bonferroni correction: Adjust the significance level
alpha <- 0.05
bonferroni_alpha <- alpha / n_tests

# Display results
#SAME - pay attention to mean severity 
bonferroni_results <- data.frame(
  Test = c("Age.Months", "Total.CBI.Score", "VineABC", "ConIN", "ConHY", 
           "CBCL.AP", "CBCL.AD.DP", "IBCL.ATT", "IBCL.AD.DP", "CSH.TSD33", 
           "SrsTotal", "SrsSCI", "SrsRRB", "SEN.TOT", "SEN.SK", "SEN.AV", 
           "SEN.SEN", "SEN.REG", "CBCL.INT", "CBCL.EXT", "mean_severity"),
  p_value = p_values,
  Bonferroni_significant = p_values < bonferroni_alpha
)

print(bonferroni_results)


#####check descriptives
#sample size
participant_count <- data %>%
  group_by(Group) %>%
  summarise(N = n())
print(participant_count)

#age
mean_age <- data %>%
  group_by(Group) %>%
  summarise(
    mean_age = mean(Age.Months, na.rm = TRUE),
    sd_age = sd(Age.Months, na.rm = TRUE)
  )
print(mean_age)


#gender
gender_count <- data%>%
  group_by(Group, Gender) %>%
  summarise(N = n())
print(gender_count)


#epilepsy
epilepsy_count <- data %>%
  group_by(Group, Epilepsy) %>%
  filter(Epilepsy == "yes") %>%  # Filter to only include those with 'yes' for epilepsy
  summarise(N = n())  # Count the number of participants
print(epilepsy_count)


#ASD diagnosis
ASD_diagnosis_count <- data %>%
  group_by(Group, ASD) %>%
  filter(ASD == "yes") %>%  # Filter to only include those with 'yes' for epilepsy
  summarise(N = n())  # Count the number of participants
print(ASD_diagnosis_count)


#ADHD diagnosis
ADHD_diagnosis_count <- data %>%
  group_by(Group, ADHD) %>%
  filter(ADHD == "yes") %>%  # Filter to only include those with 'yes' for epilepsy
  summarise(N = n())  # Count the number of participants
print(ADHD_diagnosis_count)


##check for sig differences in age/gender between case & control

#age (levenes homogeneity & t-test)
clean_age_data <- data %>%
  filter(!is.na(Age.Months))

leveneTest(Age.Months ~ Group, data = clean_age_data)

t_test_age <- t.test(Age.Months ~ Group, data = clean_age_data)
print(t_test_age)



#gender (chi square)
clean_gender_data <- data %>%
  filter(!is.na(Gender))
contingency_table <- table(data$Group, data$Gender)
print(contingency_table)

chi_square_test_gender <- chisq.test(contingency_table)
print(chi_square_test_gender)



##total CB number (chi square)
clean_CB_data <- data %>%
  filter(!is.na(Number.of.CB.reported))
contingency_table_2 <- table(clean_CB_data$Group, clean_CB_data$Number.of.CB.reported)
print(contingency_table_2)

chi_square_test_CB <- chisq.test(contingency_table_2)
print(chi_square_test_CB)

#sample size
nrow(clean_CB_data)

##try fishers exact test for CB prevalence data
fisher_test_result <- fisher.test(contingency_table_2)
print(fisher_test_result)



##CB severity by group (chi square) -- ##NOT CORRECT
clean_CB_Severity_data <- data %>%
  filter(!is.na(Total.CBI.Score))
contingency_table_3 <- table(clean_CB_data$Group, clean_CB_data$Total.CBI.Score)
print(contingency_table_3)

chi_square_test_CB_Severity <- chisq.test(contingency_table_3)
print(chi_square_test_CB_Severity)

##try fishers exact test on severity data
fisher_test_result <- fisher.test(contingency_table_3)
print(fisher_test_result)


###REDONE: checking if mean severity of CB differs by group
#using mean_severity_no_NAs --> all NA's (i.e. severity of 0) were converted to 0 for sake of statistical test 
#make sure mean_severity is created before running this - is a few steps below currently

#first - create mean_severity_no_NAs column
data <- data %>%
  mutate(mean_severity_no_NA = ifelse(is.na(mean_severity), 0, mean_severity))

#now run wilcoxon to check group differences - sig (p < 0.001)
wilcox.test(mean_severity_no_NA ~ Group, data = data)



###create column that combines IBCL & CBCL 
data <- data %>%
  mutate(CBCL.Total = coalesce(CBCL.AP, IBCL.ATT))




#######Statistics

###subset data by group (SYNGAP v TDC)
data_SYNGAP <- subset(data, Group == "SYNGAP")
data_TDC <- subset(data, Group == "TDC")

##calculate % TDC w 1 CB
participants_with_CB <- subset(
  data_TDC,
  Number.of.CB.reported >= 1
)$Participant.ID
participants_with_CB


####CREATE CBI SEVERITY AVERAGE !!! 

# Filter columns that end with 'Severity' and calculate the sum for each column
column_names <- names(data_SYNGAP)
##filter columns that end with 'severity'
severity_columns <- grep("Severity$", column_names, value = TRUE)
#average of the columns - i.e. mean IRC_Severity for each participant
mean_of_severity_columns <- data_SYNGAP %>%
  select(all_of(severity_columns)) %>%
  summarise(across(everything(), ~ mean(.x[.x > 0], na.rm = TRUE))) #change so when calculating mean, diving by number of people who are have at least a score of 1 (rather than only NAs)
print(mean_of_severity_columns)

#find ranges
range(data$Total.CBI.Score) ##CHANGE
range(data$Number.of.CB.reported)

###specifically make Total CBI (w/ average not sum)
participant_means <- data_SYNGAP %>%
  rowwise() %>%
  mutate(mean_severity = mean(c_across(all_of(severity_columns)), na.rm = TRUE)) %>%
  ungroup()

overall_CBI_severity_mean <- mean(participant_means$mean_severity, na.rm = TRUE)

print(overall_CBI_severity_mean)


##create mean severity column for all participants (both groups)
  ##adapted so mean is calculated: sum of severity scores per participant divided by number of behaviours which have a severity score of 1 or greater
data <- data %>%
  rowwise() %>%
  mutate(mean_severity = ifelse(
    all(c_across(ends_with("Severity")) <= 0 | 
          is.na(c_across(ends_with("Severity")))),
    NA,
    mean(c_across(ends_with("Severity"))[
      c_across(ends_with("Severity")) > 0
    ], na.rm = TRUE)
  )) %>%
  ungroup()

#export dataset 
#write_xlsx(data, "data_Q_CBI_with_mean_severity.xlsx")

#re-run so data_SYNGAP has mean_severity column
data_SYNGAP <- subset(data, Group == "SYNGAP")
data_TDC <- subset(data, Group == "TDC")


###START STATS

###total CB & sex
CB_sex_result <- wilcox.test(mean_severity ~ Gender, data = data_SYNGAP)
print(CB_sex_result)


###total CB & age
CB_age_result <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$Age.Months, method = "pearson")
print(CB_age_result)


##Total CB & Adaptive Behaviour 

#VineABC - significant (r = -0.53, p = 0.03)
vine_result <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$VineABC, method = "pearson")
print(vine_result)


##checking Vineland w/ SIB, PAG, STB - CHECK IF PEARSON OR SPEARMAN
vine_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$VineABC, method = "spearman")
vine_PAG <- cor.test(data_SYNGAP$PAG_Severity, data_SYNGAP$VineABC, method = "spearman")
vine_STB <- cor.test(data_SYNGAP$STB_Severity, data_SYNGAP$VineABC, method = "spearman")
print(vine_SIB)
print(vine_PAG)
print(vine_STB)

##Total CB & ADHD Traits (ConIN, ConHY, CBCL (T-test & DSM)) - n.s.

spearman_SYNGAP_ConIN <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$ConIN, method = "spearman")
spearman_SYNGAP_ConHY <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$ConHY, method = "pearson")
spearman_SYNGAP_CBCL_AP <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$CBCL.AP, method = "pearson")
spearman_SYNGAP_CBCL_AD_DP <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$CBCL.AD.DP, method = "pearson")
spearman_SYNGAP_CBCL_TOT <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$CBCL.Total, method = "pearson")

# Print results
spearman_SYNGAP_ConIN
spearman_SYNGAP_ConHY
spearman_SYNGAP_CBCL_AP ##sig (doesn't survive bonferroni corrections)
spearman_SYNGAP_CBCL_AD_DP ##sig (this one is DSM ADHD)
spearman_SYNGAP_CBCL_TOT##sig 


# Check if p-values are significant with Bonferroni correction
bonferroni_alpha <- 0.05 / 4

# Print results with adjusted significance level
cat("Spearman correlation between challenging behaviour and ConIN: p-value =", spearman_SYNGAP_ConIN$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_ConIN$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and ConHY: p-value =", spearman_SYNGAP_ConHY$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_ConHY$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and CBCL.AP: p-value =", spearman_SYNGAP_CBCL_AP$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_CBCL_AP$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and CBCL.AD.DP: p-value =", spearman_SYNGAP_CBCL_AD_DP$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_CBCL_AD_DP$p.value < bonferroni_alpha, "\n")

cat("Spearman correlation between challenging behaviour and CBCL + IBCL: p-value =", spearman_SYNGAP_CBCL_TOT$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_CBCL_TOT$p.value < bonferroni_alpha, "\n")



###checking CBCL Total score with PAG, SIB, STB
spearman_SYNGAP_CBCL_TOT_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$CBCL.Total, method = "spearman")
spearman_SYNGAP_CBCL_TOT_PAG <- cor.test(data_SYNGAP$PAG_Severity, data_SYNGAP$CBCL.Total, method = "spearman")
spearman_SYNGAP_CBCL_TOT_STB <- cor.test(data_SYNGAP$STB_Severity, data_SYNGAP$CBCL.Total, method = "spearman")

print(spearman_SYNGAP_CBCL_TOT_SIB)
print(spearman_SYNGAP_CBCL_TOT_PAG)
print(spearman_SYNGAP_CBCL_TOT_STB)


#####run test for CBCL & IBCL combined####
##(same for EC Conners?)



###Total CB & Autistic traits (SRS: Total, SCI, RRB) - all 3 sig. 
spearman_SYNGAP_SRS_T <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SrsTotal, method = "spearman")
spearman_SYNGAP_SCI <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SrsSCI, method = "spearman")
spearman_SYNGAP_RRB <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SrsRRB, method = "spearman")

# Print results
spearman_SYNGAP_SRS_T
spearman_SYNGAP_SCI
spearman_SYNGAP_RRB


# Check if p-values are significant with Bonferroni correction
bonferroni_alpha <- 0.05 / 3

# Print results with adjusted significance level
cat("Spearman correlation between challenging behaviour and SRS_T: p-value =", spearman_SYNGAP_SRS_T$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SRS_T$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and SRS SCI: p-value =", spearman_SYNGAP_SCI$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SCI$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and SRS RRB: p-value =", spearman_SYNGAP_RRB$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_RRB$p.value < bonferroni_alpha, "\n\n")


###check for SIB/PAG/STB
spearman_SRS_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SrsTotal, method = "spearman")
spearman_SRS_PAG <- cor.test(data_SYNGAP$PAG_Severity, data_SYNGAP$SrsTotal, method = "spearman")
spearman_SRS_STB <- cor.test(data_SYNGAP$STB_Severity, data_SYNGAP$SrsTotal, method = "spearman")

# Print results
spearman_SRS_SIB
spearman_SRS_PAG
spearman_SRS_STB




###Total CB & Sensory Profile scores - only sig after corrections = SEN.SK
spearman_SYNGAP_SEN_TOT <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SEN.TOT, method = "pearson")
spearman_SYNGAP_SEN_SK <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SEN.SK, method = "pearson")
spearman_SYNGAP_SEN_AV <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SEN.AV, method = "pearson")
spearman_SYNGAP_SEN_SEN <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SEN.SEN, method = "pearson")
spearman_SYNGAP_SEN_REG <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$SEN.REG, method = "spearman")

# Print results
spearman_SYNGAP_SEN_TOT #p = 0.07
spearman_SYNGAP_SEN_SK #p - 0.01
spearman_SYNGAP_SEN_AV 
spearman_SYNGAP_SEN_SEN
spearman_SYNGAP_SEN_REG


# Check if p-values are significant with Bonferroni correction
bonferroni_alpha <- 0.05 / 5

# Print results with adjusted significance level
cat("Spearman correlation between challenging behaviour and SEN Total: p-value =", spearman_SYNGAP_SEN_TOT$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_TOT$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and SEN SK: p-value =", spearman_SYNGAP_SEN_SK$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_SK$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and SEN AV: p-value =", spearman_SYNGAP_SEN_AV$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_AV$p.value < bonferroni_alpha, "\n\n")

cat("Spearman correlation between challenging behaviour and SEN SEN: p-value =", spearman_SYNGAP_SEN_SEN$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_SEN$p.value < bonferroni_alpha, "\n")

cat("Spearman correlation between challenging behaviour and SEN REG: p-value =", spearman_SYNGAP_SEN_REG$p.value, "\n")
cat("Significant with Bonferroni correction:", spearman_SYNGAP_SEN_REG$p.value < bonferroni_alpha, "\n")


###check for SIB/PAG/STB
spearman_SP_SIB <- cor.test(data_SYNGAP$SIB_Severity, data_SYNGAP$SEN.TOT, method = "spearman")
spearman_SP_PAG <- cor.test(data_SYNGAP$PAG_Severity, data_SYNGAP$SEN.TOT, method = "spearman")
spearman_SP_STB <- cor.test(data_SYNGAP$STB_Severity, data_SYNGAP$SEN.TOT, method = "spearman")

# Print results
spearman_SP_SIB
spearman_SP_PAG
spearman_SP_STB



###Total CB & Sleep -- 
result_sleep <- cor.test(data_SYNGAP$mean_severity, data_SYNGAP$CSH.TSD33, method = "pearson")
print(result_sleep)


###Total CB & CBCL Int/Ext - both sig
result_CBCL_Int <- cor.test(data$mean_severity, data$CBCL.INT, method = "pearson")
print(result_CBCL_Int)

result_CBCL_Ext <- cor.test(data$mean_severity, data$CBCL.EXT, method = "pearson")
print(result_CBCL_Ext)




#####fit REGRESSION

#remove NA's
clean_data_SYNGAP <- data_SYNGAP[complete.cases(data_SYNGAP[, c("mean_severity", "VineABC", "ConIN", "ConHY", 
                                                                "CBCL.AP", "CBCL.AD.DP", "CBCL.Total", "CSH.TSD33", 
                                                                "SrsTotal", "SrsSCI", "SrsRRB", "SEN.TOT", 
                                                                "SEN.SK", "SEN.AV", "SEN.SEN", "SEN.REG")]), ]
#create model 
regression_CB <- lm(mean_severity ~ VineABC + CBCL.AP + CSH.TSD33 + 
                      SrsTotal + SEN.TOT, data = clean_data_SYNGAP)

# Get a summary of the model
summary(regression_CB)


# Diagnostic plots
par(mfrow = c(2, 2))
plot(regression_CB)



####TRY STEPWISE REGRESSION####
# Fit a full model with all predictors
full_model <- lm(mean_severity ~ VineABC + CBCL.AP + CSH.TSD33 + SrsTotal + SEN.TOT, data = clean_data_SYNGAP)

# Apply stepwise regression (both directions)
stepwise_model <- step(full_model, direction = "backward", trace = 1)

# View the summary of the final model
summary(stepwise_model)



##PLOTS

##create custom theme
custom_theme <- theme(
  axis.title = element_text(size = 18, color = "black"),   # Axis titles
  axis.text = element_text(size = 15, color = "black"),    # Axis labels
  plot.title = element_text(size = 20, color = "black", hjust = 0.5),  # Title
  legend.title = element_text(size = 21, color = "black"), # Legend title
  legend.text = element_text(size = 19, color = "black"),   # Legend text
  plot.margin = margin(t = 10, r = 20, b = 10, l = 20), #plot margin
  panel.background = element_rect(fill = "white", color = NA),  # Set the background color to white
  panel.grid.major = element_line(color = "grey90"),  #Adjust the major grid lines color
  panel.grid.minor = element_line(color = "grey95")
)


###make bar chart for CB types 
CBI_data <- read.csv("Analysis_1.csv", header=TRUE)
head(CBI_data)

#clean data 
CBI_data$Gender <- as.factor(CBI_data$Gender)
CBI_data$Case.Control <- as.factor(CBI_data$Case.Control)

#subset data
CBI_data_SYNGAP <- subset(CBI_data, Case.Control == "SYNGAP")


#get sums for present data
present_cols <- grep("Present$", names(CBI_data_SYNGAP), value = TRUE)

# Step 2: Calculate the row-wise sum of these columns
present_sums <- colSums(CBI_data_SYNGAP[present_cols], na.rm = TRUE)

present_df <- data.frame(t(present_sums))

# Step 3: View the modified data frame
print(present_df)

###plot bar chart
#Ensure 'new_df' is in long format (if it's not already)
# Convert 'present_df' to a data frame with variables and their sums as rows
plot_data <- data.frame(
  Variable = names(present_df),
  Sum = as.numeric(present_df[1, ])
)

#change variable names
variable_names <- c("Self-Injury Behaviour", "Physical Agression", "Verbal Agression",
                    "Destruction of Property", "Anal Poking", "Stereotyped Behaviour",
                    "Inappropriate Vocalisations", "Inappropriate Removal of Clothing",
                    "Pica", "Inappropriate Sexual Behaviour", "Smearing", "Stealing", 
                    "Self-Induced Vomiting")
names(variable_names) <- names(present_df)
plot_data$Variable <- variable_names


###calculate percentages
# Calculate the total sample size
SYNGAP_sample_size <- nrow(data_SYNGAP)
print(SYNGAP_sample_size)

# Calculate the percentages and store them in a new column
plot_data <- plot_data %>%
  mutate(Percentage = (Sum / SYNGAP_sample_size) * 100)


#Plot the bar chart with variables in descending order
present_barplot <- ggplot(plot_data, aes(x = reorder(Variable, -Sum), y = Sum)) +
  geom_bar(stat = "identity", fill = "#FF7F7F") +
  labs(x = "Challenging Behaviour Type", y = "Number of Participants with CB") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 54, hjust = 1)) + 
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            vjust = -0.5, 
            color = "black",
            size = 5)  # Add percentage labels
print(present_barplot)



###make same bar chart for severity 
#get sums for severity data into new df
severity_cols <- grep("Severity$", names(CBI_data_SYNGAP), value = TRUE)
#severity_means <- colMeans(CBI_data_SYNGAP[severity_cols], na.rm = TRUE)
severity_means <- CBI_data_SYNGAP %>%
  summarise(across(all_of(severity_cols), ~ mean(.x[.x > 0], na.rm = TRUE)))
severity_df <- data.frame(severity_means)
print(severity_df)

###plot bar chart
#Ensure 'new_df' is in long format (if it's not already)
# Convert 'present_df' to a data frame with variables and their sums as rows
severity_plot_data <- data.frame(
  Variable = names(severity_df),
  Mean = as.numeric(severity_df[1, ])
)

#change variable names
variable_names <- c("Self-Injury Behaviour", "Physical Agression", "Verbal Agression",
                    "Destruction of Property", "Anal Poking", "Stereotyped Behaviour",
                    "Inappropriate Vocalisations", "Inappropriate Removal of Clothing",
                    "Pica", "Inappropriate\nSexual Behaviour", "Smearing", "Stealing", 
                    "Self-Induced Vomiting")
names(variable_names) <- names(severity_df)
severity_plot_data$Variable <- variable_names


#Plot the bar chart with variables in descending order
severity_barplot <- ggplot(severity_plot_data, aes(x = reorder(Variable, -Mean), y = Mean)) +
  geom_bar(stat = "identity", fill = "#FF7F7F") +
  labs(x = "Challenging Behaviour Type", y = "Severity of CB") +
  custom_theme +
  theme(axis.text.x = element_text(angle = 54, hjust = 1)) 
print(severity_barplot)



# Create boxplot for CB severity (means) by Group
TotalCB_Plot <- ggplot(data, aes(x = Group, y = mean_severity, fill = Group)) +
  geom_boxplot() +
  geom_point(color = "black", size = 2, width = 0.2) +
  scale_fill_manual(values = c("#FF7F7F", "#70D1C6")) +
  labs(title = "", x = "Group", y = "CB Severity Score") +
  custom_theme + 
  theme(axis.title.x = element_blank())
print(TotalCB_Plot)

##same but for CB prevalence scores
SeverityCB_Plot <- ggplot(data, aes(x = Group, y = Number.of.CB.reported, fill = Group)) +
  geom_boxplot() +
  geom_point(color = "black", size = 2, width = 0.2) +  # Add individual data points
  scale_fill_manual(values = c("#FF7F7F", "#70D1C6")) +
  labs(title = "", x = "Group", y = "Total Number of CB") +
  custom_theme + 
  theme(axis.title.x = element_blank())
print(SeverityCB_Plot)


##combine them into a single plot
combined_plot <- ggarrange(
  SeverityCB_Plot, TotalCB_Plot,
  ncol = 2, nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 16, face = "bold"),  # 👈 styling
  label.x = 0,   # left align
  label.y = 1,   # top align
  hjust = -0.5,  # tweak horizontal position
  vjust = 1.5,   # tweak vertical position
  common.legend = TRUE,
  legend = "right"
)
print(combined_plot)

##annotate to have single x axis label
combined_plot_1 <- annotate_figure(
  combined_plot,
  bottom = text_grob("Group", size = 20, color = "black", vjust = -0.2, hjust = 0.9)  # Common x-axis label
)
print(combined_plot_1)



##make a scatterplot of CB & age
Age_plot <- ggplot(data_SYNGAP, aes(x = Age.Months, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Age (months)",
       y = "Total CBI Score") +
  custom_theme
print(Age_plot)

##check if relationship is sig. - is n.s.
age_CB_correlation <- cor.test(data_SYNGAP$Age.Months, data_SYNGAP$mean_severity, method = "pearson")
print(age_CB_correlation)
#same for number of CB - n.s. 
age_noCB_correlation <- cor.test(data_SYNGAP$Age.Months, data_SYNGAP$Number.of.CB.reported, method = "spearman")
print(age_noCB_correlation)


##make a scatterplot of CB & sleep
Sleep_plot <- ggplot(data_SYNGAP, aes(x = CSH.TSD33, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "CSH Total Sleep Score",
       y = "Total CBI Score") +
  custom_theme
print(Sleep_plot)

###check if sig
sleep_CB_correlation <- cor.test(data_SYNGAP$CSH.TSD33, data_SYNGAP$mean_severity, method = "pearson")
print(sleep_CB_correlation)


#make vineABC range
vineabc_range <- range(data_SYNGAP$VineABC, na.rm = TRUE)
print(vineabc_range)

##make a scatterplot of CB & Vineland
Vine_plot <- ggplot(data_SYNGAP, aes(x = VineABC, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "",
       y = "Total CB") +
  coord_cartesian(xlim = vineabc_range) +
  custom_theme
print(Vine_plot)





###make 6 panel graph for: age, sleep, Vineland, SRS,CBCL-inattention & SP

Age_plot <- ggplot(data_SYNGAP, aes(x = Age, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Age (years)",
       y = "") +
  scale_y_continuous(limits = c(0, 30)) +
  custom_theme
print(Age_plot)

Sleep_plot <- ggplot(data_SYNGAP, aes(x = CSH.TSD33, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Child Sleep Habits Total Score",
       y = "") +
  scale_y_continuous(limits = c(0, 30)) +
  custom_theme
print(Sleep_plot)

Vine_plot <- ggplot(data_SYNGAP, aes(x = VineABC, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Vineland Adaptive Behaviour Composite",
       y = "") +
  scale_y_continuous(limits = c(0, 30)) +
  custom_theme
print(Vine_plot)


SRS_plot <- ggplot(data_SYNGAP, aes(x = SrsTotal, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Social Responsiveness Scale Total",
       y = "") +
  scale_y_continuous(limits = c(0, 30)) +
  custom_theme
print(SRS_plot)


CBCL_plot <- ggplot(data_SYNGAP, aes(x = CBCL.Total, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "CBCL Inattention T-Score",
       y = "") +
  scale_y_continuous(limits = c(0, 30)) +
  custom_theme
print(CBCL_plot)

Sensory_plot <- ggplot(data_SYNGAP, aes(x = SEN.TOT, y = mean_severity)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Sensory Profile Total Score",
       y = "") +
  scale_y_continuous(limits = c(0, 30)) +
  custom_theme
print(Sensory_plot)


##make into 6 panel plot  
combined_behaviour_plot <- ggarrange(
  Age_plot, Sleep_plot,
  Vine_plot, SRS_plot, 
  CBCL_plot, Sensory_plot, 
  ncol = 2, nrow = 3,  # Arrange in 2 columns and 2 rows
  labels = c("A", "B", "C", "D", "E", "F"),
  font.label = list(size = 16, face = "bold"),  
  label.x = 0,   # left align
  label.y = 1,   # top align
  hjust = -1.5,  # tweak horizontal position
  vjust = 1.5,   # tweak vertical position
  common.legend = TRUE,  # Common legend (if applicable)
  legend = "right"       # Place the legend on the right
)
print(combined_behaviour_plot)


# Annotate the combined plot with a common x-axis label
combined_behaviour_annotated <- annotate_figure(
  combined_behaviour_plot,
  left = text_grob("Total Challenging Behaviour Score", size = 18, color = "black", vjust = 1, hjust = 0.5, rot = 90)
)
print(combined_behaviour_annotated)


##re-create plot w/ standardised variables
# Standardize outcome (same for all plots)
data_SYNGAP$mean_severity_z <- scale(data_SYNGAP$mean_severity)

# Standardize each predictor
data_SYNGAP$Age_z        <- scale(data_SYNGAP$Age)
data_SYNGAP$Sleep_z      <- scale(data_SYNGAP$CSH.TSD33)
data_SYNGAP$Vine_z       <- scale(data_SYNGAP$VineABC)
data_SYNGAP$SRS_z        <- scale(data_SYNGAP$SrsTotal)
data_SYNGAP$CBCL_z       <- scale(data_SYNGAP$CBCL.Total)
data_SYNGAP$Sensory_z    <- scale(data_SYNGAP$SEN.TOT)


#same plots (w/ z-variables instead)
  ##this line removed from all above plots: #scale_y_continuous(limits = c(0, 30))
Age_plot_z <- ggplot(data_SYNGAP, aes(x = Age_z, y = mean_severity_z)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Age (standardized)", y = "") +
  custom_theme
Sleep_plot_z <- ggplot(data_SYNGAP, aes(x = Sleep_z, y = mean_severity_z)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Child Sleep Habits (standardized)", y = "") +
  custom_theme
Vine_plot_z <- ggplot(data_SYNGAP, aes(x = Vine_z, y = mean_severity_z)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Vineland ABC (standardized)", y = "") +
  custom_theme
SRS_plot_z <- ggplot(data_SYNGAP, aes(x = SRS_z, y = mean_severity_z)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "SRS Total (standardized)", y = "") +
  custom_theme
CBCL_plot_z <- ggplot(data_SYNGAP, aes(x = CBCL_z, y = mean_severity_z)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "CBCL Inattention (standardized)", y = "") +
  custom_theme
Sensory_plot_z <- ggplot(data_SYNGAP, aes(x = Sensory_z, y = mean_severity_z)) +
  geom_point(color = "#000080", size = 3) +
  geom_smooth(method = "lm", color = "#FF7F7F", size = 1.5) +
  labs(x = "Sensory Profile (standardized)", y = "") +
  custom_theme

#combine plots
combined_behaviour_plot <- ggarrange(
  Age_plot_z, Sleep_plot_z,
  Vine_plot_z, SRS_plot_z, 
  CBCL_plot_z, Sensory_plot_z, 
  ncol = 2, nrow = 3,  # Arrange in 2 columns and 2 rows
  labels = c("A", "B", "C", "D", "E", "F"),
  font.label = list(size = 16, face = "bold"),  
  label.x = 0,   # left align
  label.y = 1,   # top align
  hjust = -1.5,  # tweak horizontal position
  vjust = 1.5,   # tweak vertical position
  common.legend = TRUE,  # Common legend (if applicable)
  legend = "right"       # Place the legend on the right
)
print(combined_behaviour_plot)


