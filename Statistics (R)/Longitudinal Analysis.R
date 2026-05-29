library(dplyr)
library(car)
library(tidyr)
library(purrr) #for shapiro wilks
library(broom) #for shapiro wilks
library(tidyverse) 
library(lme4)
library(lmerTest)
library(ggplot2)
library(patchwork)

data <- read.csv("Longitudinal Data.csv", header=TRUE)
head(data)

#change to factors
data$Gender <- as.factor(data$Gender)
data$Group <- as.factor(data$Group)
data$ASD_1 <- as.factor(data$ASD_1)
data$ADHD_1 <- as.factor(data$ADHD_1)
data$Epilepsy_1 <- as.factor(data$Epilepsy_1)
data$ASD_2 <- as.factor(data$ASD_2)
data$ADHD_2 <- as.factor(data$ADHD_2)
data$Epilepsy_2 <- as.factor(data$Epilepsy_2)
data$ASD_3 <- as.factor(data$ASD_3)
data$ADHD_3 <- as.factor(data$ADHD_3)
data$Epilepsy_3 <- as.factor(data$Epilepsy_3)
data$Time.Point.1 <- as.factor(data$Time.Point.1)
data$Time.Point.2 <- as.factor(data$Time.Point.2)
data$Time.Point.3 <- as.factor(data$Time.Point.3)
str(data)

##remove 'space' after male
data$Gender <- trimws(data$Gender)
data$Gender <- factor(data$Gender, levels = c("Female", "Male"))
levels(data$Gender)

##same for group
data$Group <- trimws(data$Group)
data$Group <- factor(data$Group, levels = c("SYNGAP", "TDC"))
levels(data$Group)

##remove spaces & make sure each has two levels
vars_to_modify <- c("ASD_1", "ASD_2", "ASD_3", 
                    "ADHD_1", "ADHD_2", "ADHD_3",
                    "Epilepsy_1", "Epilepsy_2", "Epilepsy_3",
                    "Time.Point.1", "Time.Point.2", "Time.Point.3") # add more if needed
#loop through & make modifications
for (v in vars_to_modify) {
  
  # trim whitespace
  data[[v]] <- trimws(data[[v]])
  
  # convert empty string to NA
  data[[v]][data[[v]] == ""] <- NA
  
  # optional: standardize to "yes"/"no" lowercase
  # this does NOT change meaning, only formatting
  data[[v]] <- tolower(data[[v]])
  
  # convert to factor with allowed levels
  data[[v]] <- factor(data[[v]], levels = c("no", "yes"))
}
#check
lapply(data[vars_to_modify], levels)


##check distribution of all numeric variables 
  #RESULTS - 182 out of 337 variables are non-normally distributed; proceed with non-parametric tests generally
    #update - 204 variables non-normally distributed
#first - identify numeric columns
numeric_vars <- names(data)[sapply(data, is.numeric)]
length(numeric_vars)

#next - run shapiro wilks on all numeric columns
normality_results <- map_df(
  numeric_vars,
  \(var) {
    x <- data[[var]]
    x <- x[!is.na(x)]   # remove NA
    
    # Shapiro test requires at least 3 unique values
    if (length(unique(x)) < 3) {
      return(tibble(variable = var, p_value = NA, note = "insufficient variation"))
    }
    
    test <- shapiro.test(x)
    
    tibble(
      variable = var,
      p_value = test$p.value,
      note = ""
    )
  }
)

#adjust for multiple comparisons
normality_results <- normality_results %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr"),
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_holm = p.adjust(p_value, method = "holm")
  )

#last - identify non-normally distributed variables
non_normal <- normality_results %>%
  filter(p_fdr < 0.05)
print(non_normal, n = Inf)

#force check mean_severity columns (lots of NA's so not included above)
#BOTH NORMAL ? // lots of NAs
shapiro.test(na.omit(data$mean_severity_2))
shapiro.test(na.omit(data$mean_severity_3))


###start descriptive stats

##SAMPLE SIZE for each time point

# Specify the time point columns
timepoint_vars <- c("Time.Point.1", "Time.Point.2", "Time.Point.3")

# Make sure "yes"/"no" values are consistent (trim whitespace, lowercase)
data <- data %>%
  mutate(across(all_of(timepoint_vars), ~ tolower(trimws(.))))

# Count participants by group for each time point
sample_sizes <- data %>%
  group_by(Group) %>%
  summarise(
    T1_yes = sum(Time.Point.1 == "yes", na.rm = TRUE),
    T2_yes = sum(Time.Point.2 == "yes", na.rm = TRUE),
    T3_yes = sum(Time.Point.3 == "yes", na.rm = TRUE)
  )

sample_sizes 

#gender
# Function to get counts by group, gender, and time point
gender_by_timepoint <- map_df(timepoint_vars, function(tp) {
  data %>%
    filter(.data[[tp]] == "yes") %>%   # only include participants at this time point
    group_by(Group, Gender) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(TimePoint = tp)
})

gender_by_timepoint


###Next - age
age_vars <- c("AgeMon_1", "AgeMon_2", "AgeMon_3")

mean_age_by_group <- data %>%
  group_by(Group) %>%
  summarise(
    across(all_of(age_vars), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}")
  )

mean_age_by_group

##age in years w/ SD
age_vars_yr <- c("Age_1", "Age_2", "Age_3")

mean_age_by_group_yr <- data %>%
  group_by(Group) %>%
  summarise(
    across(
      all_of(age_vars_yr), 
      list(
        mean = ~ mean(.x, na.rm = TRUE),
         sd   = ~ sd(.x, na.rm = TRUE)
        ),
        .names = "{fn}_{col}"
      )
    )

mean_age_by_group_yr


####do the groups differ by age or gender at any time point###

#gender

#T1
t1_data <- data %>%
  filter(Time.Point.1 == "yes") %>%
  select(Group, Gender) %>%
  drop_na()
t1_table <- table(t1_data$Group, t1_data$Gender)
t1_table
fisher.test(t1_table)

#T3
t3_data <- data %>%
  filter(Time.Point.3 == "yes") %>%
  select(Group, Gender) %>%
  drop_na()
t3_table <- table(t3_data$Group, t3_data$Gender)
t3_table
fisher.test(t3_table)


###age

#T1
t1_age <- data %>%
  filter(Time.Point.1 == "yes") %>%
  select(Group, AgeMon_1) %>%
  drop_na()
wilcox.test(AgeMon_1 ~ Group, data = t1_age, exact = FALSE)

#T3
t3_age <- data %>%
  filter(Time.Point.3 == "yes") %>%
  select(Group, AgeMon_3) %>%
  drop_na()
wilcox.test(AgeMon_3 ~ Group, data = t3_age, exact = FALSE)


####START INFERENTIAL STATS####

##create custom theme
custom_theme <- theme(
  axis.title = element_text(size = 20, color = "black"),   # Axis titles
  axis.text = element_text(size = 15, color = "black"),    # Axis labels
  plot.title = element_text(size = 20, color = "black", hjust = 0.5),  # Title
  legend.title = element_text(size = 21, color = "black"), # Legend title
  legend.text = element_text(size = 19, color = "black"),   # Legend text
  plot.margin = margin(t = 10, r = 20, b = 10, l = 20), #plot margin
  panel.background = element_rect(fill = "white", color = NA),  # Set the background color to white
  panel.grid.major = element_line(color = "grey90"),  #Adjust the major grid lines color
  panel.grid.minor = element_line(color = "grey95")
)

##SRS
#first - reshape data from wide -> long format
data_long_SRS <- data %>%
  pivot_longer(
    cols = c(SrsTotal_1, SrsTotal_2, SrsTotal_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(SrsTotal))  # Keep only rows with a score

#linear mixed effects model (age x group)
model <- lmer(
  SrsTotal ~ Age * Group + (Age | Participant_ID),
  data = data_long_SRS
)

summary(model)

##plot individual trajectories
ggplot(data_long_SRS, aes(x = Age, y = SrsTotal, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme


#repeat for VineABC
data_long_VineABC <- data %>%
  pivot_longer(
    cols = c(VineABC_1, VineABC_2, VineABC_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(VineABC))   # Keep only rows with a score

#linear mixed effects model (age x group)
vine_model <- lmer(
  VineABC ~ Age * Group + (Age | Participant_ID),
  data = data_long_VineABC
)

summary(vine_model)

##plot individual trajectories (old plot)
ggplot(data_long_VineABC, aes(x = Age, y = VineABC, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  theme_minimal() +
  labs(title = "Individual VineABC Trajectories Over Age")

#new plot w/ regression lines - FIG #1
ggplot(data_long_VineABC, 
       aes(x = Age, y = VineABC, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  labs(y = "Vineland ABC Score", x = "Age (years)") +
  custom_theme


#repeat for Sensory Profile
data_long_SEN <- data %>%
  pivot_longer(
    cols = c(SP.TOT_1, SP.TOT_2, SP.TOT_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(SP.TOT))   # Keep only rows with a score

#linear mixed effects model (age x group)
SEN_model <- lmer(
  SP.TOT ~ Age * Group + (Age | Participant_ID),
  data = data_long_SEN
)

summary(SEN_model)

ggplot(data_long_SEN, aes(x = Age, y = SP.TOT, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme



##CBCL ADHD
data_combined <- data %>%
  mutate(
    CBCL_AP_combined_1 = coalesce(CBCL.AP_1, IBCL.AD.DP_1),
    CBCL_AP_combined_2 = coalesce(CBCL.AP_2, IBCL.AD.DP_2),
    CBCL_AP_combined_3 = coalesce(CBCL.AP_3, IBCL.AD.DP_3)
  )

#first reformat data into long
data_long_ADHD <- data_combined %>%
  pivot_longer(
    cols = c(CBCL_AP_combined_1, CBCL_AP_combined_2, CBCL_AP_combined_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(Time = as.integer(Time)) %>%
  filter(!is.na(CBCL_AP_combined))   # Removes rows where that timepoint doesn't exist

#linear mixed effects model (age x group)
ADHD_model <- lmer(
  CBCL_AP_combined ~ Age * Group + (Age | Participant_ID),
  data = data_long_ADHD
)

summary(ADHD_model)


ggplot(data_long_ADHD, aes(x = Age, y = CBCL_AP_combined, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme


#repeat for Sleep
data_long_sleep <- data %>%
  pivot_longer(
    cols = c(CSH.TSD33_1, CSH.TSD33_2, CSH.TSD33_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(CSH.TSD33))   # Keep only rows with a score

#linear mixed effects model (age x group)
sleep_model <- lmer(
  CSH.TSD33 ~ Age * Group + (Age | Participant_ID),
  data = data_long_sleep
)

summary(sleep_model)

ggplot(data_long_sleep, aes(x = Age, y = CSH.TSD33, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme


###CBCL- Int/Ext/Anx/Dep
##create CBCL's columns
data <- data %>%
  mutate(
    CBCL.INT.Total_1 = coalesce(CBCL.INT_1, IBCL.INT_1), #time point 1
    CBCL.EXT.Total_1 = coalesce(CBCL.EXT_1, IBCL.EXT_1),
    CBCL.ANX.Total_1 = coalesce(CBCL.AD_1, IBCL.A.D_1), 
    CBCL.INT.Total_2 = coalesce(CBCL.INT_2, IBCL.INT_2), #time point 2
    CBCL.EXT.Total_2 = coalesce(CBCL.EXT_1, IBCL.EXT_2),
    CBCL.ANX.Total_2 = coalesce(CBCL.A.D_2, IBCL.A.D_2),
    CBCL.INT.Total_3 = coalesce(CBCL.INT_3, IBCL.INT_3), #time point 3
    CBCL.EXT.Total_3 = coalesce(CBCL.EXT_3, IBCL.EXT_3),
    CBCL.ANX.Total_3 = coalesce(CBCL.A.D_3, IBCL.A.D_3)
  )


##mixed model for CBCL_INT
data_long_CBCL_INT <- data %>%
  pivot_longer(
    cols = c(CBCL.INT.Total_1, CBCL.INT.Total_2, CBCL.INT.Total_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(CBCL.INT.Total))

#model
CBCL_INT_model <- lmer(
  CBCL.INT.Total ~ Age * Group + (Age | Participant_ID),
  data = data_long_CBCL_INT
  )

summary(CBCL_INT_model)

#plot
ggplot(data_long_CBCL_INT, aes(x = Age, y = CBCL.INT.Total, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme



##mixed model for CBCL_EXT
data_long_CBCL_EXT <- data %>%
  pivot_longer(
    cols = c(CBCL.EXT.Total_1, CBCL.EXT.Total_2, CBCL.EXT.Total_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(CBCL.EXT.Total))

#model -- CHANGE 1 -> AGE??? (when all longit. participants are added)
CBCL_EXT_model <- lmer(
  CBCL.EXT.Total ~ Age * Group + (1 | Participant_ID),
  data = data_long_CBCL_EXT
)

summary(CBCL_EXT_model)

#plot
ggplot(data_long_CBCL_EXT, aes(x = Age, y = CBCL.EXT.Total, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme



##mixed model for CBCL anxiety
data_long_CBCL_ANX <- data %>%
  pivot_longer(
    cols = c(CBCL.ANX.Total_1, CBCL.ANX.Total_2, CBCL.ANX.Total_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(CBCL.ANX.Total))

#model
CBCL_ANX_model <- lmer(
  CBCL.ANX.Total ~ Age * Group + (Age | Participant_ID),
  data = data_long_CBCL_ANX
)

summary(CBCL_ANX_model)

#plot
ggplot(data_long_CBCL_ANX, aes(x = Age, y = CBCL.ANX.Total, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme




##mixed model for CBCL depression 
data_long_CBCL_DEP <- data %>%
  pivot_longer(
    cols = c(CBCL.WD_1, CBCL.W.D_2, CBCL.W.D_3,
             Age_1, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(CBCL.W.D))

#model-- CHANGE 1 -> AGE??? (when all longit. participants are added)
CBCL_DEP_model <- lmer(
  CBCL.W.D ~ Age * Group + (1 | Participant_ID),
  data = data_long_CBCL_DEP
)

summary(CBCL_DEP_model)

#plot ## FIGURE # 2
ggplot(data_long_CBCL_DEP, aes(x = Age, y = CBCL.W.D, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  labs(y = "CBCL Withdrawn/Depressed", x = "Age (years)") +
  custom_theme



##mixed model for CBI prevalence
data_long_CBI_prevalence <- data %>%
  pivot_longer(
    cols = c(CBI_Prevalence_2, CBI_Prevalence_3,
             Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(CBI_Prevalence) & !is.na(Group))  # Keep only rows with a score

data_long_CBI_prevalence_clean <- data_long_CBI_prevalence %>%
  filter(!is.na(CBI_Prevalence))

#linear mixed effects model (age x group)
CBI_prevalence_model <- lmer(
  CBI_Prevalence ~ Age * Group + (1 | Participant_ID),
  data = data_long_CBI_prevalence_clean
)

summary(CBI_prevalence_model)

ggplot(data_long_CBI_prevalence_clean, aes(x = Age, y = CBI_Prevalence, group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme


#USE MEAN SEVERITY 

##mixed model for CBI severity
data_long_mean_severity <- data %>%
  pivot_longer(
    cols = c(mean_severity_2, mean_severity_3,
             Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(!is.na(mean_severity), !is.na(Group))  # Keep only rows with a score

#linear mixed effects model (age x group)
CBI_mean_severity <- lmer(
  mean_severity ~ Age * Group + (1 | Participant_ID),
  data = data_long_mean_severity
)

summary(CBI_mean_severity)


####NEED TO ADJUST ALL NA's to ZEROs FOR THE PLOT !! 
#make separate mean_severity column where all NA's are zero's
##**note mean severity columns w/ zero don't have _ between mean and severity (naming issues)
data <- data %>%
  mutate(
    meanseverity_2_zero = replace_na(mean_severity_2, 0),
    meanseverity_3_zero = replace_na(mean_severity_3, 0)
  )

#create correct dataset (long format) with mean_severity_zero for plot
data_long_mean_severity <- data %>%
  pivot_longer(
    cols = c(meanseverity_2_zero, meanseverity_3_zero, Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(meanseverity|Age)_(\\d)(?:_zero)?"
  ) %>%
  mutate(Time = as.integer(Time))

#plot
ggplot(data_long_mean_severity,
       aes(x = Age, y = meanseverity,
           group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    aes(group = Group, color = Group),
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    alpha = 0.25
  ) +
  scale_color_discrete(na.translate = FALSE) +
  custom_theme

#need to remove NA's from plot


##mixed model for CBI_SIB 
#convert data to long WITH FILTERING OUT 0's IN SIB
data_long_SIB_Severity <- data %>%
  pivot_longer(
    cols = c(SIB_Severity_2, SIB_Severity_3,
             Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(
    !is.na(SIB_Severity),
    SIB_Severity > 0
  )
#convert group to factor
data_long_SIB_Severity <- data_long_SIB_Severity %>%
  mutate(
    Group = factor(Group),
    Participant_ID = factor(Participant_ID)
  )
#filter for only SYNGAP group (as no TDC scored on this)
data_long_SIB_Severity_SYNGAP <- data_long_SIB_Severity %>%
  filter(Group == "SYNGAP") %>%
  mutate(Participant_ID = droplevels(Participant_ID))


#run model
SIB_Severity_model <- lmer(
  SIB_Severity ~ Age + (1 | Participant_ID),
  data = data_long_SIB_Severity_SYNGAP
)

summary(SIB_Severity_model)



#UPDATED PLOT (to handle only 1 group) - FIG #3!!! 
ggplot(data_long_SIB_Severity_SYNGAP,
       aes(x = Age, y = SIB_Severity,
           group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(group = 1),
    linewidth = 0.8,
    alpha = 0.25
  ) +
  labs(y = "CBI Self-Injury Score", x = "Age (years)") +
  custom_theme +
  guides(color = "none")


##mixed model for CBI_PAG
#convert data to long WITH FILTERING OUT 0's IN PAG
data_long_PAG_Severity <- data %>%
  pivot_longer(
    cols = c(PAG_Severity_2, PAG_Severity_3,
             Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(
    !is.na(PAG_Severity),
    PAG_Severity > 0
  )
#convert group to factor
data_long_PAG_Severity <- data_long_PAG_Severity %>%
  mutate(
    Group = factor(Group),
    Participant_ID = factor(Participant_ID)
  )
#filter for only SYNGAP group (as no TDC scored on this)
data_long_PAG_Severity_SYNGAP <- data_long_PAG_Severity %>%
  filter(Group == "SYNGAP") %>%
  mutate(Participant_ID = droplevels(Participant_ID))


#run model
PAG_Severity_model <- lmer(
  PAG_Severity ~ Age + (1 | Participant_ID),
  data = data_long_PAG_Severity_SYNGAP
)

summary(PAG_Severity_model)

#UPDATED PLOT (to handle only 1 group)
ggplot(data_long_PAG_Severity_SYNGAP,
       aes(x = Age, y = PAG_Severity,
           group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(group = 1),
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme +
  guides(color = "none")

##mixed model for CBI_STB
#convert data to long WITH FILTERING OUT 0's IN STB
data_long_STB_Severity <- data %>%
  pivot_longer(
    cols = c(STB_Severity_2, STB_Severity_3,
             Age_2, Age_3),
    names_to = c(".value", "Time"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(
    Time = as.integer(Time)
  ) %>%
  filter(
    !is.na(STB_Severity),
    STB_Severity > 0
  )
#convert group to factor
data_long_STB_Severity <- data_long_STB_Severity %>%
  mutate(
    Group = factor(Group),
    Participant_ID = factor(Participant_ID)
  )
#filter for only SYNGAP group (as no TDC scored on this)
data_long_STB_Severity_SYNGAP <- data_long_STB_Severity %>%
  filter(Group == "SYNGAP") %>%
  mutate(Participant_ID = droplevels(Participant_ID))


#run model
STB_Severity_model <- lmer(
  STB_Severity ~ Age + (1 | Participant_ID),
  data = data_long_STB_Severity_SYNGAP
)

summary(STB_Severity_model)

#UPDATED PLOT (to handle only 1 group)
ggplot(data_long_STB_Severity_SYNGAP,
       aes(x = Age, y = STB_Severity,
           group = Participant_ID, color = Group)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 1.0) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(group = 1),
    linewidth = 0.8,
    alpha = 0.25
  ) +
  custom_theme +
  guides(color = "none")




##make trajectory plot - FIG 2 !!! 
####multi-panel plot for all n.s. findings 

make_trajectory_plot <- function(data, y_var, y_label) {
  ggplot(
    data,
    aes(
      x = Age,
      y = .data[[y_var]],
      group = Participant_ID,
      color = Group
    )
  ) +
    geom_line(alpha = 0.8) +
    geom_point(alpha = 1.0) +
    geom_smooth(
      aes(group = Group, color = Group),
      method = "lm",
      se = TRUE,
      linewidth = 0.8,
      alpha = 0.25
    ) +
    scale_color_discrete(na.translate = FALSE) +
    labs(x = NULL, y = y_label) +
    custom_theme
}

##make each plot - here SRS
p_SRS <- make_trajectory_plot(
  data_long_SRS,
  "SrsTotal",
  "Autistic Traits"
)

#CBCL inattention
p_ADHD <- make_trajectory_plot(
  data_long_ADHD,
  "CBCL_AP_combined",
  "CBCL\nInattention"
)

#CBCL internalising
p_CBCL_INT <- make_trajectory_plot(
  data_long_CBCL_INT,
  "CBCL.INT.Total",
  "Internalising"
)

#CBCL externalising
p_CBCL_EXT <- make_trajectory_plot(
  data_long_CBCL_EXT,
  "CBCL.EXT.Total",
  "Externalising"
)

#CBCL anxiety
p_CBCL_ANX <- make_trajectory_plot(
  data_long_CBCL_ANX,
  "CBCL.ANX.Total",
  "Anxiety"
)

#CBCL depression
p_CBCL_DEP <- make_trajectory_plot(
  data_long_CBCL_DEP,
  "CBCL.W.D",
  "Depression"
)

#sensory profile
p_Sensory <- make_trajectory_plot(
  data_long_SEN,
  "SP.TOT",
  "Sensory Profile"
)

#Sleep
p_Sleep <- make_trajectory_plot(
  data_long_sleep,
  "CSH.TSD33",
  "Sleep Habits"
)

#CBI prevalence
p_CBI_prev <- make_trajectory_plot(
  data_long_CBI_prevalence_clean,
  "CBI_Prevalence",
  "CBI Prevalence"
)

#CBI severity
p_CBI_sev <- make_trajectory_plot(
  data_long_mean_severity,
  "meanseverity",
  "CBI Severity"
)

#create 2x4 panel plot 
multi_panel_plot <-
  (p_SRS | p_ADHD) / 
  (p_Sensory | p_Sleep) /
  (p_CBI_prev | p_CBI_sev)
#display
multi_panel_plot

#add A-F labels & shared X-axis label
final_plot <-
  multi_panel_plot +
  plot_layout(guides = "collect") +
  theme(
    legend.position = "right",          
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.95)
  ) +
  plot_annotation(
    tag_levels = "A",
    caption = "Age (years)"
  ) &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.17, 1.0),
    plot.caption = element_text(
      size = 18,
      hjust = 0.5,
      margin = margin(t = 10)
    )
  )

final_plot

#removed plots
#(p_CBCL_INT | p_CBCL_EXT) /
#  (p_CBCL_ANX | p_CBCL_DEP) / 


#add panel labels
multi_panel_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.95),
    legend.position = "right"
  )

