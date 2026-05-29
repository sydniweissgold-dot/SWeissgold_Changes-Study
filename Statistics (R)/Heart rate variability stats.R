library(dplyr)
library(car)
library(tidyr)
library(ARTool)
library(ggplot2)

data <- read.csv("0000_CHANGES_HRV_OUTPUT.csv", header=TRUE)

head(data)

#check & change to factors
is.factor((data$Gender))
is.factor((data$Group))
data$Gender <- as.factor(data$Gender)
data$Group <- as.factor(data$Group)

#change to numeric
# List of columns to convert
cols_to_numeric <- c("Age", "baseline.Recording_Length..seconds.", "baseline.Number_of_Beats_Rejected",
                     "arousing.Recording_Length..seconds.", "arousing.Number_of_Beats_Rejected",
                     "calming.Recording_Length..seconds.", "calming.Number_of_Beats_Rejected")

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
numeric_data <- data[, sapply(data, is.numeric)][, -1]
normality_tests <- sapply(numeric_data, function(x) shapiro.test(x)$p.value)

alpha <- 0.05  #significance level
bonferroni_threshold <- alpha / length(normality_tests)  # Adjusted alpha

normality_results <- data.frame(Variable = names(normality_tests),
                                P_Value = normality_tests,
                                Normal_Distribution = normality_tests > bonferroni_threshold)

print(normality_results)



##check for sig differences in age/gender between case & control

#age 
leveneTest(Age ~ Group, data = data) 

t_test_age <- t.test(Age ~ Group, data = data)
print(t_test_age) 


#gender 
contingency_table <- table(data$Group, data$Gender)
print(contingency_table)

chi_square_test_gender <- chisq.test(contingency_table)
print(chi_square_test_gender)


##basic t-tests
t_test <- t.test(baseline.bpm ~ Group, data = data)
print(t_test)



####convert data to long & remove NA's

data_long <- data %>%
  pivot_longer(
    cols = -c(Participant.ID, Group, Gender, Age),  # Keep participant/group variables
    names_to = c("Condition", "Metric"),  # Split variable names into two parts
    names_sep = "\\.",  # Assumes variable names are formatted like "rmssd_baseline"
    values_to = "Value"  # Stores the actual data values
  )
head(data_long)


data_cleaned <- data_long %>%
  group_by(Participant.ID) %>%
  filter(!any(is.na(Value))) %>%  # Remove participants with any NA in any metric
  ungroup()  # Ungroup to prevent issues in later analyses

head(data_cleaned)



###check whether we need to control for age, beats rejected, recording length

#age
cor_age_results <- data_cleaned %>%
  group_by(Metric) %>%
  summarise(correlation = cor.test(Age, Value, method = "spearman", use = "pairwise.complete.obs")$estimate,
            p_value = cor.test(Age, Value, method = "spearman", use = "pairwise.complete.obs")$p.value)

print(cor_age_results, n = 23)


##beats rejected
#reshape data for beats rejected
beats_rejected_long <- data %>%
  select(Participant.ID, starts_with("baseline.Number_of_Beats_Rejected"),
         starts_with("arousing.Number_of_Beats_Rejected"),
         starts_with("calming.Number_of_Beats_Rejected")) %>%
  pivot_longer(
    cols = -Participant.ID,
    names_to = "Condition",
    names_pattern = "(.*)\\.Number_of_Beats_Rejected",
    values_to = "Beats_Rejected"
  )


#merge with cleaned data
data_with_rejected <- data_cleaned %>%
  left_join(beats_rejected_long, by = c("Participant.ID", "Condition"))


#run spearman w beats rejected and all metrics
cor_rejected_results <- data_with_rejected %>%
  group_by(Condition, Metric) %>%
  summarise(
    correlation = cor.test(Beats_Rejected, Value, method = "spearman", use = "pairwise.complete.obs")$estimate,
    p_value = cor.test(Beats_Rejected, Value, method = "spearman", use = "pairwise.complete.obs")$p.value,
    .groups = "drop"  # Optional: prevents nested grouping in output
  )

print(cor_rejected_results, n = Inf)




##now check recording length

#reshape recording length
recording_length_long <- data %>%
  select(Participant.ID, starts_with("baseline.Recording_Length"),
         starts_with("arousing.Recording_Length"),
         starts_with("calming.Recording_Length")) %>%
  pivot_longer(
    cols = -Participant.ID,
    names_to = "Condition",
    names_pattern = "(.*)\\.Recording_Length.*",  # Adjust pattern if needed
    values_to = "Recording_Length"
  )


#merge with long-format data
data_with_length <- data_with_rejected %>%  # assuming you already have beats rejected joined
  left_join(recording_length_long, by = c("Participant.ID", "Condition"))

#run spearman correlations
cor_length_results <- data_with_length %>%
  group_by(Condition, Metric) %>%
  summarise(
    correlation = cor.test(Recording_Length, Value, method = "spearman", use = "pairwise.complete.obs")$estimate,
    p_value = cor.test(Recording_Length, Value, method = "spearman", use = "pairwise.complete.obs")$p.value,
    .groups = "drop"
  )

print(cor_length_results, n = Inf)



###Multivariate Stats###

data_cleaned$Participant.ID <- as.factor(data_cleaned$Participant.ID)
data_cleaned$Condition <- as.factor(data_cleaned$Condition)
data_cleaned$Group <- as.factor(data_cleaned$Group)



##Run Art ANOVA for specific list of variables

selected_metrics <- c("bpm", "ibi", "rmssd", "sdnn", "pnn50", "vlf", 
                      "lf", "hf", "s", "sd1", "sd2") 

anova_results <- lapply(selected_metrics, function(m) {
  model <- art(Value ~ Group * Condition + (1 | Participant.ID), 
               data = data_cleaned %>% filter(Metric == m))
  
  return(list(Metric = m, ANOVA = anova(model)))
})

anova_results



###run Art ANOVA for variables individually w/ plots

##bpm
art_model_bpm <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "bpm"))
anova(art_model_bpm)

data_bpm <- data_cleaned %>%
  filter(Metric == "bpm")
ggplot(data = data_bpm, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "bpm", title = "bpm by Group and Condition")


###RMSSD
art_model_rmssd <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "rmssd"))
anova(art_model_rmssd)

data_rmssd <- data_cleaned %>%
  filter(Metric == "rmssd")
ggplot(data = data_rmssd, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "rmssd", title = "rmssd by Group and Condition")

#get means
means_rmssd <- data_rmssd %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_rmssd)


##SDNN
art_model_sdnn <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "sdnn"))
anova(art_model_sdnn)

data_sdnn <- data_cleaned %>%
  filter(Metric == "sdnn")
ggplot(data = data_sdnn, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "sdnn", title = "sdnn by Group and Condition")

#get means
means_sdnn <- data_sdnn %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_sdnn)


#pNN50
art_model_pnn50 <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "pnn50"))
anova(art_model_pnn50)

data_pnn50 <- data_cleaned %>%
  filter(Metric == "pnn50")
ggplot(data = data_pnn50, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "pnn50", title = "pnn50 by Group and Condition")

#get means
means_pnn50 <- data_pnn50 %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_pnn50)

#VLF
art_model_vlf <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "vlf"))
anova(art_model_vlf)

data_vlf <- data_cleaned %>%
  filter(Metric == "vlf")
ggplot(data = data_vlf, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "vlf", title = "vlf by Group and Condition")


#LF
art_model_lf <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "lf"))
anova(art_model_lf)

data_lf <- data_cleaned %>%
  filter(Metric == "lf")
ggplot(data = data_lf, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "lf", title = "lf by Group and Condition")


#HF
art_model_hf <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "hf"))
anova(art_model_hf)

data_hf <- data_cleaned %>%
  filter(Metric == "hf")
ggplot(data = data_hf, aes(x = Condition, y = Value, fill = Group)) + geom_boxplot() + labs(y = "hf", title = "hf by Group and Condition")

#get means:
means_hf <- data_hf %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_hf)

#S
art_model_s <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "s"))
anova(art_model_s)

#SD1
art_model_sd1 <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "sd1"))
anova(art_model_sd1)

#SD2
art_model_sd2 <- art(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned %>% filter(Metric == "sd2"))
anova(art_model_sd2)




















