library(dplyr)
library(car)
library(tidyr)
library(ez) #for ez anova
library(rstatix) #for mixed anova
library(ggplot2)

data <- read.csv("0000_CHANGES_PUPIL_OUTPUT.csv", header=TRUE)

#checks on data
head(data)

#check & change to factors
is.factor((data$Gender))
is.factor((data$Group))
data$Gender <- as.factor(data$Gender)
data$Group <- as.factor(data$Group)
data$Best.Eye <- as.factor(data$Best.Eye)

#change to numeric (most were int)
# List of columns to convert
cols_to_numeric <- c("Age", "Best.Eye.Valid.Samples.Pre", "Best.Eye.Valid.Samples.Post",
                     "Total.Blinks", "Baseline.Blinks", "Arousing.Blinks", "Calming.Blinks", 
                     "Sample.Rate..Hz.", "Derivative.Threshold..mm.s.")

# Convert columns to numeric
data[cols_to_numeric] <- lapply(data[cols_to_numeric], function(x) as.numeric(as.character(x)))

# Verify the changes
str(data)  # Check the structure to confirm successful conversion


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
    mean_age = mean(Age, na.rm = TRUE),
    sd_age = sd(Age, na.rm = TRUE)
  )
print(mean_age)


#gender
gender_count <- data%>%
  group_by(Group, Gender) %>%
  summarise(N = n())
print(gender_count)


##Check for normal distribution##
###check distribution & do bonferroni correction for multiple comparisons

#results - ###variables that need NON-parametric tests:
  #change from arousing to calming
  #calming blinks & total blinks
  #all duration (in min)
  #pre & post butterworth metrics


# Select only numeric columns (excluding the first column)
numeric_data <- data[, sapply(data, is.numeric)][, -1]

# Run Shapiro-Wilk test only on valid columns (with >3 non-NA unique values)
normality_tests <- sapply(numeric_data, function(x) {
  x <- na.omit(x)
  if (length(unique(x)) > 3) {
    return(shapiro.test(x)$p.value)
  } else {
    return(NA)  # skip columns with too few or identical values
  }
})

#see which columns were skipped
  #was just i.e. sample rate - analysis metrics that were skipped
skipped_vars <- names(numeric_data)[sapply(numeric_data, function(x) length(unique(na.omit(x))) <= 3)]
cat("Skipped columns (identical or insufficient variation):", skipped_vars, "\n")

# Ensure it's a data frame with matching length
normality_results <- data.frame(
  Variable = names(normality_tests),
  Shapiro_p = as.numeric(normality_tests)
)

# Bonferroni correction
alpha <- 0.05  #significance level
bonferroni_threshold <- alpha / length(normality_tests)  # Adjusted alpha

# Output results
normality_results <- data.frame(Variable = names(normality_tests),
                                P_Value = normality_tests,
                                Normal_Distribution = normality_tests > bonferroni_threshold)

print(normality_results)


##check for sig differences in age/gender between case & control

#age (levenes homogeneity & t-test)
leveneTest(Age ~ Group, data = data) #ns

t_test_age <- t.test(Age ~ Group, data = data)
print(t_test_age) #ns


#gender (chi square)
contingency_table <- table(data$Group, data$Gender)
print(contingency_table)

chi_square_test_gender <- chisq.test(contingency_table)
print(chi_square_test_gender) #n.s.


###check if there are sig differences in recording length between groups
  #all n.s.
# List of your three variables
metrics <- c("Baseline.Duration..min.", 
             "Arousing.Duration..min.", 
             "Calming.Duration..min.")

# Run Wilcoxon rank-sum tests for each variable
wilcoxon_results <- lapply(metrics, function(var) {
  test <- wilcox.test(data[[var]] ~ data$Group, exact = FALSE)
  data.frame(
    Variable = var,
    W = test$statistic,
    p_value = test$p.value
  )
})

# Combine into one data frame
wilcoxon_results <- do.call(rbind, wilcoxon_results)

print(wilcoxon_results)


###now - check to see if pupil size is affected by age
  #all n.s. 
# List of pupil metrics
pupil_metrics <- c("Baseline.Mean", "Arousing.Mean", "Calming.Mean")

# Compute Pearson correlation for each
cor_results <- lapply(pupil_metrics, function(var) {
  test <- cor.test(data[[var]], data$Age, method = "pearson")
  data.frame(
    Variable = var,
    Correlation = test$estimate,
    p_value = test$p.value
  )
})

# Combine into one data frame
cor_results <- do.call(rbind, cor_results)

# Apply Bonferroni correction for 3 tests
cor_results$Bonferroni_p <- p.adjust(cor_results$p_value, method = "bonferroni")

# Add significance flag
cor_results$Significance <- ifelse(cor_results$Bonferroni_p < 0.05, "Significant", "Not Significant")

print(cor_results)


##next - see if number of blinks differs by condition &/or group
  #results - blink count differs by condition (but not group)
    #further investigations - baseline differs from calm / arousing groups (but blinks during calm / arousing don't differ from each other)

#first re-shape data into long format
# Gather the blink measures
data_long <- data %>%
  select(Participant.ID, Group, Baseline.Blinks, Arousing.Blinks, Calming.Blinks) %>%
  pivot_longer(
    cols = c(Baseline.Blinks, Arousing.Blinks, Calming.Blinks),
    names_to = "Condition",
    values_to = "Blinks"
  )

# Clean up condition names
data_long$Condition <- factor(data_long$Condition,
                              levels = c("Baseline.Blinks", "Arousing.Blinks", "Calming.Blinks"),
                              labels = c("Baseline", "Arousing", "Calming"))
#anova (this is parametric - can handle mild non-normality - i.e. from calming blinks variable)
anova_results <- ezANOVA(
  data = data_long,
  dv = .(Blinks),
  wid = .(Participant.ID),
  within = .(Condition),
  between = .(Group),
  type = 3,
  detailed = TRUE
)

print(anova_results)

##re-run blink rate anova - excluding baseline (for reporting in thesis)
# Filter out baseline
data_blinks_no_baseline <- data_long %>%
  filter(Condition %in% c("Calming", "Arousing"))

# Run ANOVA
anova_no_baseline <- ezANOVA(
  data = data_blinks_no_baseline,
  dv = .(Blinks),
  wid = .(Participant.ID),
  within = .(Condition),
  between = .(Group),
  type = 3,
  detailed = TRUE
)

print(anova_no_baseline)


#check for outliers in pupil size variables - found none! 

# List of pupil size variables
pupil_metrics <- c("Baseline.Mean", "Arousing.Mean", "Calming.Mean")

# Function to flag outliers (+/- 3 SD)
find_outliers <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  lower <- mean_x - 3*sd_x
  upper <- mean_x + 3*sd_x
  return(which(x < lower | x > upper))  # return row indices
}

# Apply to 3 metrics
outlier_indices <- lapply(pupil_metrics, function(var) {
  find_outliers(data[[var]])
})

names(outlier_indices) <- pupil_metrics

# Print outliers
outlier_indices


##PRIMARY STATS - if there is a difference in pupil size by group (or condition)
#n.s.

#first - reshape data into long format
data_long_2 <- data %>%
  select(Participant.ID, Group, Baseline.Mean, Arousing.Mean, Calming.Mean) %>%
  pivot_longer(
    cols = c(Baseline.Mean, Arousing.Mean, Calming.Mean),
    names_to = "Condition",
    values_to = "PupilSize"
  )

# Clean up condition names
data_long_2$Condition <- factor(data_long_2$Condition,
                              levels = c("Baseline.Mean", "Arousing.Mean", "Calming.Mean"),
                              labels = c("Baseline", "Arousing", "Calming"))

#run mixed-anova
anova_pupil_results <- data_long_2 %>%
  anova_test(dv = PupilSize, wid = Participant.ID, within = Condition, between = Group)

# View ANOVA table
get_anova_table(anova_pupil_results)


##re-run - only looking at differences between arousing & calm (not baseline)
  #n.s.
# Keep only Arousing and Calming
data_subset <- data %>%
  select(Participant.ID, Group, Arousing.Mean, Calming.Mean) %>%
  pivot_longer(
    cols = c(Arousing.Mean, Calming.Mean),
    names_to = "Condition",
    values_to = "PupilSize"
  )

# Clean up condition names
data_subset$Condition <- factor(data_subset$Condition,
                                levels = c("Arousing.Mean", "Calming.Mean"),
                                labels = c("Arousing", "Calming"))

#run anova
anova_pupil_2_results <- data_subset %>%
  anova_test(dv = PupilSize, wid = Participant.ID, within = Condition, between = Group)

# View results
get_anova_table(anova_pupil_2_results)


#plot this
ggplot(data_subset, aes(x = Condition, y = PupilSize, color = Group, group = Group)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(0.1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, position = position_dodge(0.1)) +
  stat_summary(fun = mean, geom = "point", size = 3, position = position_dodge(0.1)) +
  theme_minimal() +
  labs(title = "",
       y = "Mean Pupil Size", x = "Condition")



#see if there a difference in the change between calming / arousing conditions between groups
  #n.s.
t_test_diff <- t.test(data$`Arousing...Calming` ~ data$Group)
t_test_diff

