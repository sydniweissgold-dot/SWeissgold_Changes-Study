library(dplyr)
library(car)
library(tidyr)
library(Hmisc) #for correlation matrices


data_hrv <- read.csv("0000_CHANGES_HRV_OUTPUT_MANUAL.csv", header=TRUE)
head(data_hrv)

data_pupil <- read.csv("0000_CHANGES_PUPIL_OUTPUT.csv", header=TRUE)
head(data_pupil)

#remove underscore in pupil particpant ID (to match hrv)
data_pupil$Participant.ID <- gsub("_", "", data_pupil$Participant.ID)

##merge datasets by Participant ID (keep age, gender, group once in new dataset)
data_combined <- data_hrv %>%
  full_join(data_pupil, by = "Participant.ID") %>%
  mutate(
    Age = coalesce(Age.x, Age.y),
    Gender = coalesce(Gender.x, Gender.y),
    Group = coalesce(Group.x, Group.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y"))

head(data_combined)

#re-order
data_combined <- data_combined %>%
  select(Participant.ID, Age, Gender, Group, everything())


data_combined$Gender <- as.factor(data_combined$Gender)
data_combined$Group <- as.factor(data_combined$Group)
data_combined$Best.Eye <- as.factor(data_combined$Best.Eye)

cols_to_numeric <- c("Age", "baseline.NN50", "calming.NN50", "arousing.NN50", 
                     "Lux", "Best.Eye.Valid.Samples.Pre", "Best.Eye.Valid.Samples.Post",
                     "Total.Blinks", "Baseline.Blinks", "Arousing.Blinks", "Calming.Blinks", 
                     "Sample.Rate..Hz.", "Derivative.Threshold..mm.s.")
data_combined[cols_to_numeric] <- lapply(data_combined[cols_to_numeric], function(x) as.numeric(as.character(x)))
str(data_combined) 


#first - filter out participants wiht NA for either baselinermssd/arousingmean
sample_sizes <- data_combined %>%
  filter(!is.na(baseline.rmssd) & !is.na(Arousing.Mean)) %>%
  group_by(Group) %>%
  summarise(n = n())
#calculate sample size (by group)
sample_sizes


#get descriptives for included group
filtered_data <- data_combined %>%
  filter(!is.na(baseline.rmssd) & !is.na(Arousing.Mean))

age_summary <- filtered_data %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    age_mean = mean(Age, na.rm = TRUE),
    age_sd = sd(Age, na.rm = TRUE),
    age_min = min(Age, na.rm = TRUE),
    age_max = max(Age, na.rm = TRUE)
  )

gender_summary <- filtered_data %>%
  group_by(Group, Gender) %>%
  summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = Gender, values_from = n, values_fill = 0)


age_summary
gender_summary




###variables that need non-parametric tests (taken from hrv/pupil stats code)
  #hrv: mostly non-parametric 
  #pupil: lux, arousing from calming, total blinks, calming blinks, 
          #all durations


#run correlation for all baseline variables (mean pupil size & hrv)
vars_to_correlate <- c("baseline.ibi", "baseline.bpm", "baseline.rmssd", 
                       "baseline.sdnn", "baseline.pnn50", "baseline.sd1", 
                       "baseline.sd2", "baseline.vlf", "baseline.lf", "baseline.hf")

physio_data <- data_combined %>%
  select(Baseline.Mean, all_of(vars_to_correlate))

physio_matrix <- as.matrix(physio_data)

corr_results <- rcorr(physio_matrix, type = "spearman")

corr_long <- data.frame(
  Variable1 = "baseline.mean",
  Variable2 = vars_to_correlate,
  Correlation = corr_results$r["Baseline.Mean", vars_to_correlate],
  P_value = corr_results$P["Baseline.Mean", vars_to_correlate]
)

corr_long



#same for arousing
arousing_vars <- c("arousing.ibi", "arousing.bpm", "arousing.rmssd", 
                   "arousing.sdnn", "arousing.pnn50", "arousing.sd1", 
                   "arousing.sd2", "arousing.vlf", "arousing.lf", "arousing.hf")

arousing_data <- data_combined %>%
  select(Arousing.Mean, all_of(arousing_vars))

arousing_matrix <- as.matrix(arousing_data)

arousing_corr <- rcorr(arousing_matrix, type = "spearman")

arousing_corr_long <- data.frame(
  Variable1 = "Arousing.Mean",
  Variable2 = arousing_vars,
  Correlation = arousing_corr$r["Arousing.Mean", arousing_vars],
  P_value = arousing_corr$P["Arousing.Mean", arousing_vars]
)

arousing_corr_long



##same for calming
calming_vars <- c("calming.ibi", "calming.bpm", "calming.rmssd", 
                  "calming.sdnn", "calming.pnn50", "calming.sd1", 
                  "calming.sd2", "calming.vlf", "calming.lf", "calming.hf")

calming_data <- data_combined %>%
  select(Calming.Mean, all_of(calming_vars))

calming_matrix <- as.matrix(calming_data)

calming_corr <- rcorr(calming_matrix, type = "spearman")

calming_corr_long <- data.frame(
  Variable1 = "Calming.Mean",
  Variable2 = calming_vars,
  Correlation = calming_corr$r["Calming.Mean", calming_vars],
  P_value = calming_corr$P["Calming.Mean", calming_vars]
)

calming_corr_long



