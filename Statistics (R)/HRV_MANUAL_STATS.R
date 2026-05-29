library(dplyr)
library(car)
library(tidyr)
library(ARTool)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(stringr)
library(lme4)
library(lmerTest)
library(rstatix)

data <- read.csv("0000_CHANGES_HRV_OUTPUT_MANUAL.csv", header=TRUE)

#checks on data
head(data)

#check & change to factors
is.factor((data$Gender))
is.factor((data$Group))
data$Gender <- as.factor(data$Gender)
data$Group <- as.factor(data$Group)

#change to numeric
# List of columns to convert
cols_to_numeric <- c("Age")


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

#results - mixed TRUE / FALSE - will assume non-parametric but need to remember to check! 

# Select only numeric columns (excluding the first column)
numeric_data <- data[, sapply(data, is.numeric)][, -1]
# Run Shapiro-Wilk test on all numeric columns (apart from first - participant ID)
normality_tests <- sapply(numeric_data, function(x) shapiro.test(x)$p.value)

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



##look at age x BPM
#####need to run outlier filter first!! 
#filter data
data_bpm_arousing <- data_cleaned_no_outliers %>%
  filter(Metric == "bpm", Condition == "arousing")

#run linear regression
lm_bpm_arousing <- lm(Value ~ Age * Group, data = data_bpm_arousing)
summary(lm_bpm_arousing)

#plot data
ggplot(data_bpm_arousing, aes(x = Age, y = Value, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "BPM (Arousing) by Age and Group",
       x = "Age",
       y = "BPM") +
  theme_minimal()



#gender (chi square)
contingency_table <- table(data$Group, data$Gender)
print(contingency_table)

chi_square_test_gender <- chisq.test(contingency_table)
print(chi_square_test_gender) #n.s.


##basic t-tests
t_test <- t.test(baseline.bpm ~ Group, data = data)
print(t_test)



####convert data to long & remove NA's (clean data)
#necessary for multivariate stats



# Step 1: Clean column names
names(data) <- names(data) %>%
  str_replace_all("\\.", "_") %>%       # replace all dots with underscores
  str_replace_all("_+", "_") %>%        # collapse multiple underscores
  str_replace_all("_$", "")             # remove trailing underscores



# Step 2: Pivot longer
data_long <- data %>%
  pivot_longer(
    cols = -c(Participant_ID, Group, Gender, Age),
    names_to = c("Condition", "Metric"),
    names_pattern = "([^_]*)_(.*)",   # split at first underscore
    values_to = "Value"
  )


# Check the structure
head(data_long)

# Remove participants with ANY missing values across all dependent variables
data_cleaned <- data_long %>%
  group_by(Participant_ID) %>%
  filter(!any(is.na(Value))) %>%  # Remove participants with any NA in any metric
  ungroup()  # Ungroup to prevent issues in later analyses


# Check structure
head(data_cleaned)



###REMOVE OUTLIERS +/- 3SD from mean###
remove_outliers <- function(df) {
  df %>%
    group_by(Group, Condition, Metric) %>%  # remove outliers within participant
    mutate(
      mean_val = mean(Value, na.rm = TRUE),
      sd_val = sd(Value, na.rm = TRUE),
      is_outlier = abs(Value - mean_val) > 3 * sd_val,
      is_outlier = ifelse(is.na(is_outlier), FALSE, is_outlier)  # keep rows where SD couldn't be computed
    ) %>%
    filter(!is_outlier) %>%
    select(-mean_val, -sd_val, -is_outlier) %>%
    ungroup()
}


# Apply outlier removal
data_cleaned_no_outliers <- remove_outliers(data_cleaned)

#print how many data points were removed
cat("Removed", nrow(data_cleaned) - nrow(data_cleaned_no_outliers), "outliers.\n")


# Count how many participants were removed bc of outliers
participants_with_outliers <- setdiff(unique(data_cleaned$Participant_ID), 
                                      unique(data_cleaned_no_outliers$Participant_ID))

cat("Participants with *all* their data removed due to outliers:", length(participants_with_outliers), "\n")
print(participants_with_outliers)



##number of participants who had at least one dataset removed from outliers
affected_participants <- data_cleaned %>%
  anti_join(data_cleaned_no_outliers, by = c("Participant_ID", "Condition", "Metric", "Value")) %>%
  distinct(Participant_ID)

cat("Participants with at least one outlier removed:", nrow(affected_participants), "\n")
#list participant IDs
print(affected_participants)


###SET CONDITION VARIABLE ORDER###
#define consistent condition ordering
conditions_to_keep <- c("baseline", "calming", "arousing")

# Ensure Condition always has the same order
data_cleaned_no_outliers$Condition <- factor(
  data_cleaned_no_outliers$Condition,
  levels = conditions_to_keep
)

# Create a "no-baseline" version for plots/analyses
data_no_baseline <- data_cleaned_no_outliers %>%
  filter(Condition %in% c("calming", "arousing")) %>%
  droplevels()


####SET JITTER FOR PLOTS LATER
  #if want jitter = box_jitter
  #if want boxplot (no jitter) = geom_boxplot()
box_jitter <- list(
  geom_boxplot(outlier.shape = NA, alpha = 0.6),
  geom_jitter(
    aes(color = Group),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    alpha = 0.7,
    size = 2
  )
)


###Multivariate Stats###


##make sure group & condition are factors
data_cleaned_no_outliers$Participant_ID <- as.factor(data_cleaned_no_outliers$Participant_ID)
data_cleaned_no_outliers$Condition <- as.factor(data_cleaned_no_outliers$Condition)
data_cleaned_no_outliers$Group <- as.factor(data_cleaned_no_outliers$Group)


##ART MODEL -- for non-parametric & repeated measures 
##investigate - what does the '1' mean in art model 

##Run Art ANOVA for specific list of variables

# Define the specific list of metrics you want to analyze
selected_metrics <- c("bpm", "ibi", "rmssd", "sdnn", "pnn50", "vlf", 
                      "lf", "hf", "sd1", "sd2") 


#test to see if random effect is needed: 
#results - show lots of variance in rmssd between participants so random effect is suitable
#library(lme4)
#test_model <- lmer(Value ~ Group * Condition + (1 | Participant.ID), data = data_cleaned_no_outliers %>% filter(Metric == "rmssd"))
#summary(test_model)



###RUN ANOVAs w/OUT BASELINE
# Shared filter to exclude 'baseline'
conditions_to_keep <- c("arousing", "calming")


## bpm
data_bpm <- data_no_baseline %>%
  filter(Metric == "bpm")

art_model_bpm <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_bpm)
anova(art_model_bpm)

plot_bpm_no_baseline <- ggplot(data = data_bpm, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming",
               "arousing" = "Arousing")
  ) +
  labs(y = "BPM", x = NULL)
plot_bpm_no_baseline

###bpm plot w/ added jitter
plot_bpm_no_baseline_jitter <- ggplot(data = data_bpm, aes(x = Condition, y = Value, fill = Group)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # hide boxplot outlier dots so they don't duplicate
  geom_jitter(aes(color = Group),   # color points by group
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              alpha = 0.7, size = 2) +
  labs(y = "bpm", title = "BPM") + 
  theme_minimal()
plot_bpm_no_baseline_jitter

#get means
means_bpm <- data_bpm %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_bpm)


## ibi
data_ibi <- data_no_baseline %>%
  filter(Metric == "ibi")

art_model_ibi <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_ibi)
anova(art_model_ibi)

plot_ibi_no_baseline <- ggplot(data = data_ibi, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming",
               "arousing" = "Arousing")
  ) +
  labs(y = "IBI", x = NULL)
plot_ibi_no_baseline

#get means
means_ibi <- data_ibi %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_ibi)


## rmssd
data_rmssd <- data_no_baseline %>%
  filter(Metric == "rmssd")

art_model_rmssd_no_baseline <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_rmssd)
anova(art_model_rmssd_no_baseline)

plot_rmssd_no_baseline <- ggplot(data = data_rmssd, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "RMSSD", x = NULL)
plot_rmssd_no_baseline

means_rmssd <- data_rmssd %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_rmssd)


## sdnn
data_sdnn <- data_no_baseline %>%
  filter(Metric == "sdnn")

art_model_sdnn <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_sdnn)
anova(art_model_sdnn)

plot_sdnn_no_baseline <- ggplot(data = data_sdnn, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "SDNN", x = NULL)
plot_sdnn_no_baseline

#get means
means_sdnn <- data_sdnn %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_sdnn)


## pnn50
data_pnn50 <- data_no_baseline %>%
  filter(Metric == "pnn50")

art_model_pnn50 <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_pnn50)
anova(art_model_pnn50)

plot_pnn50_no_baseline <- ggplot(data = data_pnn50, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "pNN50", x = NULL)
plot_pnn50_no_baseline

means_pnn50 <- data_pnn50 %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_pnn50)


## vlf
data_vlf <- data_no_baseline %>%
  filter(Metric == "vlf")

art_model_vlf <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_vlf)
anova(art_model_vlf)

plot_vlf_no_baseline <- ggplot(data = data_vlf, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "VLF", x = "")
plot_vlf_no_baseline


#get means
means_vlf <- data_vlf %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_vlf)


## lf
data_lf <- data_no_baseline %>%
  filter(Metric == "lf")

art_model_lf <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_lf)
anova(art_model_lf)

plot_lf_no_baseline <- ggplot(data = data_lf, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "LF", x = "")
plot_lf_no_baseline


#get means
means_lf <- data_lf %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')

print(means_lf)


## hf
data_hf <- data_no_baseline %>%
  filter(Metric == "hf")

art_model_hf <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_hf)
anova(art_model_hf)

plot_hf_no_baseline <- ggplot(data = data_hf, aes(x = Condition, y = Value, fill = Group)) + 
  box_jitter + 
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "HF", x = "")
plot_hf_no_baseline

#get means
means_hf <- data_hf %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_hf)


## log HF
data_log_hf <- data_hf %>%
  mutate(HFLogValue = log(Value))

art_model_log_hf <- art(HFLogValue ~ Group * Condition + (1 | Participant_ID), data = data_log_hf)
anova(art_model_log_hf)

####parametric anova - linear mixed effects model (as anova not great w/ repeated measures)
parametric_log_hf <- lmer(HFLogValue ~ Group * Condition + (1 | Participant_ID), 
                     data = data_log_hf)
anova(parametric_log_hf)



means_log_hf <- data_log_hf %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(HFLogValue, na.rm = TRUE), .groups = 'drop')
print(means_log_hf)

plot_log_hf_no_baseline <- ggplot(data_log_hf, aes(x = Condition, y = HFLogValue, fill = Group)) +
  box_jitter +
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "LogHF", x = "")
plot_log_hf_no_baseline


## sd1
data_sd1 <- data_cleaned_no_outliers %>%
  filter(Metric == "sd1", Condition %in% conditions_to_keep)

art_model_sd1 <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_sd1)
anova(art_model_sd1)

means_sd1 <- data_sd1 %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_sd1)

plot_sd1_no_baseline <- ggplot(data_sd1, aes(x = Condition, y = Value, fill = Group)) +
  box_jitter +
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "SD1")
plot_sd1_no_baseline


## sd2
data_sd2 <- data_cleaned_no_outliers %>%
  filter(Metric == "sd2", Condition %in% conditions_to_keep)

art_model_sd2 <- art(Value ~ Group * Condition + (1 | Participant_ID), data = data_sd2)
anova(art_model_sd2)

means_sd2 <- data_sd2 %>%
  group_by(Group, Condition) %>%
  summarize(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
print(means_sd2)

plot_sd2_no_baseline <- ggplot(data_sd2, aes(x = Condition, y = Value, fill = Group)) +
  box_jitter +
  scale_x_discrete(
    labels = c("calming" = "Calming", #capitalise calming and arousing
               "arousing" = "Arousing")
  ) +
  labs(y = "SD2")
plot_sd2_no_baseline


###make 9-panel plot
#plot_bpm_no_baseline <- plot_bpm_no_baseline
#plot_ibi_no_baseline <- plot_ibi_no_baseline
#plot_rmssd_no_baseline <- plot_rmssd_no_baseline
#plot_sdnn_no_baseline <- plot_sdnn_no_baseline
#plot_pnn50_no_baseline <- plot_pnn50_no_baseline
#plot_vlf_no_baseline <- plot_vlf_no_baseline
#plot_lf_no_baseline <- plot_lf_no_baseline
#plot_hf_no_baseline <- plot_hf_no_baseline
#plot_log_hf_no_baseline <- plot_log_hf_no_baseline
#plot_sd1_no_baseline <- plot_sd1_no_baseline
#plot_sd2_no_baseline <- plot_sd2_no_baseline

# Combine plots with shared legend 
#combined_plot_no_baseline <-
# (plot_bpm_no_baseline | plot_ibi_no_baseline | plot_rmssd_no_baseline) / #row 1
# (plot_sdnn_no_baseline | plot_pnn50_no_baseline | plot_vlf_no_baseline) / #row 2
# (plot_lf_no_baseline | plot_hf_no_baseline | plot_log_hf_no_baseline) / #row 3
# (plot_spacer() | plot_sd1_no_baseline | plot_sd2_no_baseline | plot_spacer()) + #row 4 (spacer to centre)
# plot_layout(guides = "collect") &
# theme(legend.position = "right") &
# labs(x = "Condition")

# Show the plot
#combined_plot_no_baseline



####LOG TRANSFORM HRV METRICS
log_metrics <- c(
  "bpm", "ibi", "sdnn", "pnn50",
  "rmssd", "lf", "vlf", "hf",
  "sd1", "sd2"
)

# Create empty lists to store outputs
anova_results <- list()
model_results <- list()
means_results <- list()
plot_results <- list()


# Loop through each metric
for(metric_name in log_metrics){
  
  # Filter metric
  metric_data <- data_no_baseline %>%
    filter(Metric == metric_name) %>%
    
    # Log transform
    mutate(
      LogValue = log(Value)
    )
  
  # Linear mixed model
  #model <- lmer(
    #LogValue ~ Group * Condition + (1 | Participant_ID),
    #data = metric_data
  #)
  
  # Repeated-measures ANOVA
  model <- anova_test(
    data = metric_data,
    dv = LogValue,
    wid = Participant_ID,
    within = Condition,
    between = Group
  )
  
  # Store model
  model_results[[metric_name]] <- model
  
  # ANOVA table
  anova_table <- get_anova_table(model)
  
  anova_results[[metric_name]] <- anova_table
  
  # Print results
  cat("\n============================\n")
  cat("Metric:", metric_name, "\n")
  
  print(anova_table)
  
  # Means
  means <- metric_data %>%
    group_by(Group, Condition) %>%
    summarise(
      mean_value = mean(LogValue, na.rm = TRUE),
      sd_value = sd(LogValue, na.rm = TRUE),
      .groups = "drop"
    )
  
  means_results[[metric_name]] <- means
  
  print(means)
}



##data checks e.g. RMSSD here
data_RMSSD_check <- data_no_baseline %>%
  filter(Metric == "rmssd") %>%
  mutate(LogValue = log(Value + 0.001))

# Raw values
ggplot(data_RMSSD_check, aes(x = Value)) +
  geom_histogram(bins = 30) +
  labs(title = "Raw RMSSD")

# Log-transformed values
ggplot(data_RMSSD_check, aes(x = LogValue)) +
  geom_histogram(bins = 30) +
  labs(title = "Log RMSSD")



###run ART ANOVA for log transformed pNN50 - NOT USING (NOT SURE IF ITS CORRECT)
# Filter pNN50 data
data_pNN50 <- data_no_baseline %>%
  dplyr::filter(Metric == "pnn50") %>%
  mutate(
    LogValue = log(Value),
    Participant_ID = factor(Participant_ID),
    Group = factor(Group),
    Condition = factor(Condition)
  )
# ART ANOVA on log-transformed values
model_pNN50_art <- art(
  LogValue ~ Group * Condition + (1 | Participant_ID),
  data = data_pNN50
)
# ANOVA table
anova_pNN50_art <- anova(model_pNN50_art)
print(anova_pNN50_art)




##specify results for each metric (interchange) -example here is hf
anova_results[["hf"]]
summary(model_results[["hf"]])
plot_results[["hf"]]
means_results[["hf"]]



#### MAKE PLOTS FOR EACH METRIC

for(metric_name in log_metrics){
  
  # Filter metric and log transform
  metric_data <- data_no_baseline %>%
    filter(Metric == metric_name) %>%
    mutate(
      LogValue = log(Value)
    )
  
  # Create plot
  p <- ggplot(
    metric_data,
    aes(
      x = Condition,
      y = LogValue,
      fill = Group
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.6
    ) +
    
    geom_jitter(
      aes(color = Group),
      position = position_jitterdodge(
        jitter.width = 0.2,
        dodge.width = 0.75
      ),
      alpha = 0.7,
      size = 2
    ) +
    
    scale_x_discrete(
      labels = c(
        "calming" = "Calming",
        "arousing" = "Arousing"
      )
    ) +
    
    labs(
      y = paste("Ln", toupper(metric_name)),
      x = NULL
    ) +
    
    theme_minimal()
  
  # Store plot
  plot_results[[metric_name]] <- p
  
  # Print plot
  print(p)
}



#HR plot
plot_results[["bpm"]] + plot_results[["ibi"]]

#time domain plot
plot_results[["rmssd"]] + plot_results[["sdnn"]] + plot_results[["pnn50"]]

#frequency domain plot
plot_results[["vlf"]] + plot_results[["lf"]] + plot_results[["hf"]]

#non-linear plot
plot_results[["sd1"]] + plot_results[["sd2"]]




#### HISTOGRAMS OF LOG-TRANSFORMED METRICS

# Create empty list for histograms
histogram_results <- list()

for(metric_name in log_metrics){
  
  # Filter metric and log transform
  metric_data <- data_no_baseline %>%
    filter(Metric == metric_name) %>%
    mutate(
      LogValue = log(Value)
    )
  
  # Histogram
  p_hist <- ggplot(metric_data, aes(x = LogValue)) +
    
    geom_histogram(
      bins = 30,
      color = "black",
      fill = "skyblue",
      alpha = 0.7
    ) +
    
    labs(
      title = paste("Histogram of Log", toupper(metric_name)),
      x = paste("Log", toupper(metric_name)),
      y = "Frequency"
    ) +
    
    theme_minimal()
  
  # Store histogram
  histogram_results[[metric_name]] <- p_hist
  
  # Print histogram
  print(p_hist)
}



##Make a Poincare plot

# Read the file
rr_data <- read.csv("14CHA_cleaned_peaks.arousing.csv", header = FALSE)

# Assign column names (adjust if your actual file has headers)
colnames(rr_data) <- c("PeakIndex", "Timestamp", "RR")
# Ensure RR is numeric
rr_data$RR <- as.numeric(rr_data$RR)

# Check structure
head(rr_data)


# Create RR_n and RR_n+1
poincare_data <- rr_data %>%
  mutate(RR_n = lag(RR),
         RR_n1 = RR) %>%
  filter(!is.na(RR_n))  # Remove first row with NA


# Create the Poincaré plot
ggplot(poincare_data, aes(x = RR_n, y = RR_n1)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  labs(x = expression(RR[n]~"(ms)"),
       y = expression(RR[n+1]~"(ms)")) +
  theme_minimal()




##making arousing & calming plots for 1 participant: 14CHA
# Load and label each condition
load_rr_data <- function(file_path, condition_label) {
  df <- read.csv(file_path, header = FALSE)
  colnames(df) <- c("PeakIndex", "Timestamp", "RR")
  df$RR <- as.numeric(df$RR)
  df$Condition <- condition_label
  return(df)
}

# Load both files
rr_arousing <- load_rr_data("14CHA_cleaned_peaks.arousing.csv", "Arousing")
rr_calming <- load_rr_data("14CHA_cleaned_peaks.calming.csv", "Calming")

# Function to create Poincaré plot
make_poincare_plot <- function(data, title_text) {
  poincare_data <- data %>%
    mutate(RR_n = lag(RR),
           RR_n1 = RR) %>%
    filter(!is.na(RR_n))
  
  ggplot(poincare_data, aes(x = RR_n, y = RR_n1)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "steelblue4") +
    labs(title = title_text,
         x = expression(RR[n]~"(ms)"),
         y = expression(RR[n+1]~"(ms)")) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5, margin = margin(t = 3), size = 14, face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)
    )
}

# Create both plots
plot_calming <- make_poincare_plot(rr_calming, "Calming")
plot_arousing <- make_poincare_plot(rr_arousing, "Arousing")


# Combine side-by-side
(plot_calming + plot_arousing) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.02, 0.95))




####MAKE PLOTS (w/out baseline)



###MAKE BPM & IBI PLOT####


#PLOT
combined_HR_plot <- ggarrange(
  plot_bpm_no_baseline,
  plot_ibi_no_baseline,
  ncol = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B") 
)

combined_HR_plot <- annotate_figure(
  combined_HR_plot,
  bottom = text_grob("Condition", face = "bold", size = 12)
)

combined_HR_plot


###MAKE TIME-DOMAIN PLOT - w/out baseline####

#make sure condition labels are capitalised
plot_rmssd_no_baseline <- plot_rmssd_no_baseline

plot_sdnn_no_baseline <- plot_sdnn_no_baseline

plot_pnn50_no_baseline <- plot_pnn50_no_baseline

#PLOT
combined_time_domain_plot <- ggarrange(
  plot_rmssd_no_baseline,
  plot_sdnn_no_baseline,
  plot_pnn50_no_baseline,
  ncol = 3,
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B", "C")
)

combined_time_domain_plot <- annotate_figure(
  combined_time_domain_plot,
  bottom = text_grob("Condition", face = "bold", size = 12)
)

combined_time_domain_plot




###MAKE FREQUENCY-DOMAIN PLOT####

#make sure condition labels are capitalised
plot_vlf_no_baseline <- plot_vlf_no_baseline

plot_lf_no_baseline <- plot_lf_no_baseline

plot_hf_no_baseline <- plot_hf_no_baseline

plot_log_hf_no_baseline <- plot_log_hf_no_baseline



#PLOT
combined_frequency_domain_plot <- ggarrange(
  plot_vlf_no_baseline,
  plot_lf_no_baseline,
  plot_hf_no_baseline,
  plot_log_hf_no_baseline,
  ncol = 2, nrow = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B", "C", "D")
)

combined_frequency_domain_plot <- annotate_figure(
  combined_frequency_domain_plot,
  bottom = text_grob("Condition", face = "bold", size = 12)
)

combined_frequency_domain_plot




#####plots - trade so group is on X-axis and condition is different blocks
##NOT USING
#bpm
plot_bpm_no_baseline_jitter <- ggplot(data = data_bpm, 
                                      aes(x = Group, y = Value, fill = Condition)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6,
               position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = Condition),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.7, size = 2) +
  labs(y = "BPM", title = "BPM", x = NULL) + 
  theme_minimal() +
  scale_fill_discrete(
    labels = c("calming" = "Calming",
               "arousing" = "Arousing")
  ) +
  scale_color_discrete(
    labels = c("calming" = "Calming",
               "arousing" = "Arousing")
  )

plot_bpm_no_baseline_jitter

