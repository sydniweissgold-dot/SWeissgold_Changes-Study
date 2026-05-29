library(dplyr)
library(car)
library(tidyr)
library(purrr) #for shapiro wilks
library(broom) #for shapiro wilks
library(tidyverse) 
library(lme4)
library(lmerTest)
library(ggplot2)
library(stringr)
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
  
  data[[v]] <- trimws(data[[v]])
  
  data[[v]][data[[v]] == ""] <- NA
  
  data[[v]] <- tolower(data[[v]])
  
  data[[v]] <- factor(data[[v]], levels = c("no", "yes"))
}
lapply(data[vars_to_modify], levels)


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



##combining CBCL elements (if needed?)
data <- data %>%
  mutate(
    CBCL_AP_combined_1 = coalesce(CBCL.AP_1, IBCL.AD.DP_1),
    CBCL_AP_combined_2 = coalesce(CBCL.AP_2, IBCL.AD.DP_2),
    CBCL_AP_combined_3 = coalesce(CBCL.AP_3, IBCL.AD.DP_3)
  )




####PERSISTENCE ANALYSIS

#starting w/ SIB T2 & T3 (SYNGAP only)
#filter for SYNGAP only
data_SYNGAP <- data %>%
  filter(Group == "SYNGAP")

#create 3 groups
data_SYNGAP <- data_SYNGAP %>%
  mutate(
    SIB_persistence = case_when(
      SIB_Prevalence_2 == 0 & SIB_Prevalence_3 == 0 ~ "Absent",
      SIB_Prevalence_2 == 1 & SIB_Prevalence_3 == 1 ~ "Persistent",
      SIB_Prevalence_2 != SIB_Prevalence_3            ~ "Transient",
      TRUE ~ NA_character_
    ),
    SIB_persistence = factor(SIB_persistence,
                             levels = c("Absent", "Transient", "Persistent"))
  )
#check results
SIB_summary <- data_SYNGAP %>%
  filter(!is.na(SIB_persistence)) %>%   # remove NA's (all people who are missing T3 data)
  count(SIB_persistence) %>%
  mutate(
    percent = round(100 * n / sum(n), 1)
  )
SIB_summary

#plot of prevalence categories
ggplot(SIB_summary, aes(x = SIB_persistence, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(percent, "%")), vjust = -0.5) +
  custom_theme +
  labs(y = "Count", x = "SIB Persistence")

##if wanting to see who NA's were
data_SYNGAP %>%
  filter(is.na(SIB_persistence)) %>%
  select(Participant_ID, SIB_Prevalence_2, SIB_Prevalence_3)


#SAME for PAG
data_SYNGAP <- data_SYNGAP %>%
  mutate(
    PAG_persistence = case_when(
      PAG_Prevalence_2 == 0 & PAG_Prevalence_3 == 0 ~ "Absent",
      PAG_Prevalence_2 == 1 & PAG_Prevalence_3 == 1 ~ "Persistent",
      PAG_Prevalence_2 != PAG_Prevalence_3            ~ "Transient",
      TRUE ~ NA_character_
    ),
    PAG_persistence = factor(PAG_persistence,
                             levels = c("Absent", "Transient", "Persistent"))
  )
#check results 
PAG_summary <- data_SYNGAP %>%
  filter(!is.na(PAG_persistence)) %>%   # remove NA's (all people who are missing T3 data)
  count(PAG_persistence) %>%
  mutate(
    percent = round(100 * n / sum(n), 1)
  )
PAG_summary
#plot of prevalence categories
ggplot(PAG_summary, aes(x = PAG_persistence, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(percent, "%")), vjust = -0.5) +
  custom_theme +
  labs(y = "Count", x = "PAG Persistence")


#SAME for STB 
data_SYNGAP <- data_SYNGAP %>%
  mutate(
    STB_persistence = case_when(
      STB_Prevalence_2 == 0 & STB_Prevalence_3 == 0 ~ "Absent",
      STB_Prevalence_2 == 1 & STB_Prevalence_3 == 1 ~ "Persistent",
      STB_Prevalence_2 != STB_Prevalence_3            ~ "Transient",
      TRUE ~ NA_character_
    ),
    STB_persistence = factor(STB_persistence,
                             levels = c("Absent", "Transient", "Persistent"))
  )
#check results 
STB_summary <- data_SYNGAP %>%
  filter(!is.na(STB_persistence)) %>%   # remove NA's (all people who are missing T3 data)
  count(STB_persistence) %>%
  mutate(
    percent = round(100 * n / sum(n), 1)
  )
STB_summary
#plot of prevalence categories
ggplot(STB_summary, aes(x = STB_persistence, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(percent, "%")), vjust = -0.5) +
  custom_theme +
  labs(y = "Count", x = "STB Persistence")


###try for CBCL Item 18

##first - collapse 1 & 2's
data_SYNGAP <- data_SYNGAP %>%
  mutate(
    CBCL_SIB_1_coll = ifelse(CBCL_SIB_1 >= 1, 1,
                             ifelse(CBCL_SIB_1 == 0, 0, NA)),
    CBCL_SIB_2_coll = ifelse(CBCL_SIB_2 >= 1, 1,
                             ifelse(CBCL_SIB_2 == 0, 0, NA)),
    CBCL_SIB_3_coll = ifelse(CBCL_SIB_3 >= 1, 1,
                             ifelse(CBCL_SIB_3 == 0, 0, NA))
  )

##next - only looking at T2 & T3 (to check reliability w/ CBI)
data_SYNGAP <- data_SYNGAP %>%
  mutate(
    CBCL_SIB_persistence = case_when(
      CBCL_SIB_2_coll == 0 & CBCL_SIB_3_coll == 0 ~ "absent",
      CBCL_SIB_2_coll == 1 & CBCL_SIB_3_coll == 1 ~ "persistent",
      CBCL_SIB_2_coll != CBCL_SIB_3_coll ~ "transient",
      TRUE ~ NA_character_
    ),
    CBCL_SIB_persistence = factor(CBCL_SIB_persistence,
                             levels = c("absent", "transient", "persistent"))
  )
#check results 
CBCL_SIB_summary_2T <- data_SYNGAP %>%
  filter(!is.na(CBCL_SIB_persistence)) %>%   # remove NA's (all people who are missing T3 data)
  count(CBCL_SIB_persistence) %>%
  mutate(
    percent = round(100 * n / sum(n), 1)
  )
CBCL_SIB_summary_2T
#plot of prevalence categories
ggplot(CBCL_SIB_summary_2T, aes(x = CBCL_SIB_persistence, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(percent, "%")), vjust = -0.5) +
  custom_theme +
  labs(y = "Count", x = "CBCL_SIB Persistence 2T")


#different logic as have 3 time points now 
  ##only including participants in this if they have 2 or more timepoints
data_SYNGAP <- data_SYNGAP %>%
  mutate(
    # number of observed timepoints
    CBCL_SIB_n_obs = rowSums(!is.na(cbind(CBCL_SIB_1_coll,
                                          CBCL_SIB_2_coll,
                                          CBCL_SIB_3_coll))),
    # number of endorsed timepoints
    CBCL_SIB_n_yes = rowSums(cbind(CBCL_SIB_1_coll,
                                   CBCL_SIB_2_coll,
                                   CBCL_SIB_3_coll) == 1,
                             na.rm = TRUE),
    
    CBCL_SIB_persistence = case_when(
      CBCL_SIB_n_obs >= 2 & CBCL_SIB_n_yes == 0 ~ "absent",
      CBCL_SIB_n_obs >= 2 & CBCL_SIB_n_yes == CBCL_SIB_n_obs ~ "persistent",
      CBCL_SIB_n_obs >= 2 & CBCL_SIB_n_yes > 0 & CBCL_SIB_n_yes < CBCL_SIB_n_obs ~ "transient",
      TRUE ~ NA_character_
    ),
    
    CBCL_SIB_persistence = factor(
      CBCL_SIB_persistence,
      levels = c("absent", "transient", "persistent")
    )
  )

#check results
CBCL_SIB_summary <- data_SYNGAP %>%
  filter(!is.na(CBCL_SIB_persistence)) %>%
  count(CBCL_SIB_persistence) %>%
  mutate(
    percent = round(100 * n / sum(n), 1)
  )

CBCL_SIB_summary

#plot 
ggplot(CBCL_SIB_summary, aes(x = CBCL_SIB_persistence, y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(percent, "%")), vjust = -0.5) +
  custom_theme +
  labs(y = "Count", x = "CBCL_SIB 3T Persistence")



###Kruskall-wallis - 
  #relationships between persistence of SIB w/ other factors

##Factors to examine 
variables <- c(
  "Age_2",
  "VineABC_2", "VineCOMM_2", "VineDLS_2", "VineSOC_2",
  "SrsRRB_2","SrsSCI_2",
  "CBCL_AP_combined_2",
  "CSH.TSD33_2",
  "SP.REG_2", "SP.AV_2", "SP.SEN_2", "SP.SK_2"
)


#run test for SIB_Persistence
SIB_results <- lapply(variables, function(v) {
  formula <- as.formula(paste(v, "~ SIB_persistence"))
  test <- kruskal.test(formula, data = data_SYNGAP)
  list(variable = v, test = test)
})

SIB_results


#check direction of VineABC - FIG 4
ggplot(data_SYNGAP %>% filter(!is.na(SIB_persistence)),
       aes(x = SIB_persistence, y = VineABC_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "SIB Persistence",
       y = "Vineland ABC Score") +
  custom_theme



#run for PAG_Persistence
PAG_results <- lapply(variables, function(v) {
  formula <- as.formula(paste(v, "~ PAG_persistence"))
  test <- kruskal.test(formula, data = data_SYNGAP)
  list(variable = v, test = test)
})

PAG_results

##run for STB
STB_results <- lapply(variables, function(v) {
  formula <- as.formula(paste(v, "~ STB_persistence"))
  test <- kruskal.test(formula, data = data_SYNGAP)
  list(variable = v, test = test)
})

STB_results

##check direction of STB sig results: Vine ABC, COMM, SOC, SP REG, SEN, SK

#check direction of VineABC
p_VineABC <- ggplot(data_SYNGAP %>% filter(!is.na(STB_persistence)),
       aes(x = STB_persistence, y = VineABC_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "",
       y = "Vineland ABC") +
  custom_theme


#check direction of VineCOMM
p_VineCOMM <- ggplot(data_SYNGAP %>% filter(!is.na(STB_persistence)),
       aes(x = STB_persistence, y = VineCOMM_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "",
       y = "Vineland\nCommunication") +
  custom_theme


#check direction of VineSOC
p_VineSOC <- ggplot(data_SYNGAP %>% filter(!is.na(STB_persistence)),
       aes(x = STB_persistence, y = VineSOC_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "",
       y = "Vineland\nSocialisation") +
  custom_theme


#check direction of SP REG
p_SP_REG <- ggplot(data_SYNGAP %>% filter(!is.na(STB_persistence)),
       aes(x = STB_persistence, y = SP.REG_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "",
       y = "Sensory Prorfile\nRegistation") +
  custom_theme


#check direction of SP SEN
p_SP_SEN <- ggplot(data_SYNGAP %>% filter(!is.na(STB_persistence)),
       aes(x = STB_persistence, y = SP.SEN_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "",
       y = "Sensory Profile\nSensitivity") +
  custom_theme



#check direction of SP SK 
p_SP_SK <- ggplot(data_SYNGAP %>% filter(!is.na(STB_persistence)),
       aes(x = STB_persistence, y = SP.SK_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, fill = "lightcoral") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7, color = "lightcoral") +
  labs(x = "",
       y = "Sensory Profile\nSeeking") +
  custom_theme



##make STB multi-panel plot -- FIG 5 !!! 

#first - set plot coordinates
p_VineABC  <- p_VineABC  + coord_cartesian(ylim = c(20, 80))
p_VineCOMM <- p_VineCOMM + coord_cartesian(ylim = c(20, 80))
p_VineSOC  <- p_VineSOC  + coord_cartesian(ylim = c(20, 80))

p_SP_REG <- p_SP_REG + coord_cartesian(ylim = c(20, 80))
p_SP_SEN <- p_SP_SEN + coord_cartesian(ylim = c(20, 80))
p_SP_SK  <- p_SP_SK  + coord_cartesian(ylim = c(20, 80))

#combine plots into one figure
STB_multi_panel_plot <-
  (p_VineABC | p_SP_REG) /
  (p_VineCOMM | p_SP_SEN) /
  (p_VineSOC | p_SP_SK)

STB_multi_panel_plot


STB_multi_panel_plot +
  plot_annotation(tag_levels = "A", 
  caption = "Stereotyped Behaviour Persistence") &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.17, 1.0),
    plot.caption = element_text(
      size = 18,      # ← change font size here
      hjust = 0.5     # center the shared x-axis label
    )
  )


###not currently reporting this 
###example binary logistic regression - using only persistent & absent (transient not included?)

#SIB (DV) & age as predictor 
   ##--> adjust age as predictor based on kruskall-wallis findings
#first - make binary data set for absent & persistent only (no transient)
data_SIB_binary <- data_SYNGAP %>%
  filter(SIB_persistence %in% c("Persistent", "Absent")) %>%
  mutate(
    SIB_persistent_binary = ifelse(SIB_persistence == "Persistent", 1, 0)
  )
#run model
glm(SIB_persistent_binary ~ VineABC_2, 
    data = data_SIB_binary, family = binomial) |> summary()



####MULTIVARIATE ANALYSES across ALL behav. factors

###model doesn't seem to be a good fit 

#SIB first
glm(SIB_persistent_binary ~ VineABC_2 + SrsTotal_2,
    data = data_SIB_binary, family = binomial) |> summary()

#lm w/ vif for multicollinearity
vif(
  lm(
    VineABC_2 ~ SP.TOT_2 + SrsTotal_2 + CBCL_AP_combined_2 + CSH.TSD33_2,
    data = data_SIB_binary
  )
)


##same for PAG 
#make PAG binary data
data_PAG_binary <- data_SYNGAP %>%
  filter(PAG_persistence %in% c("persistent", "absent")) %>%
  mutate(
    PAG_persistent_binary = ifelse(PAG_persistence == "persistent", 1, 0)
  )
#run regression
glm(PAG_persistent_binary ~ VineABC_2 + SP.TOT_2 + SrsTotal_2 + 
      CBCL_AP_combined_2 + CSH.TSD33_2,
    data = data_PAG_binary, family = binomial) |> summary()

#lm w/ vif for multicollinearity
#vif > 5 indicates multicollinearity 
vif(
  lm(
    VineABC_2 ~ SP.TOT_2 + SrsTotal_2 + CBCL_AP_combined_2 + CSH.TSD33_2,
    data = data_PAG_binary
  )
)


#SAME FOR STB 
data_STB_binary <- data_SYNGAP %>%
  filter(STB_persistence %in% c("persistent", "absent")) %>%
  mutate(
    STB_persistent_binary = ifelse(STB_persistence == "persistent", 1, 0)
  )
#run regression
glm(STB_persistent_binary ~ VineABC_2 + SP.TOT_2 + SrsTotal_2 + 
      CBCL_AP_combined_2 + CSH.TSD33_2,
    data = data_STB_binary, family = binomial) |> summary()

#lm w/ vif for multicollinearity
#vif > 5 indicates multicollinearity 
vif(
  lm(
    VineABC_2 ~ SP.TOT_2 + SrsTotal_2 + CBCL_AP_combined_2 + CSH.TSD33_2,
    data = data_STB_binary
  )
)



##ADD ANS MEASURES

#set-up
data_hrv <- read.csv("0000_CHANGES_HRV_OUTPUT_MANUAL.csv", header=TRUE)
head(data_hrv)

data_pupil <- read.csv("0000_CHANGES_PUPIL_OUTPUT.csv", header=TRUE)
head(data_pupil)

data_pupil$Participant.ID <- gsub("_", "", data_pupil$Participant.ID)

data$Participant_ID <- data$Participant_ID |>
  str_extract("^\\d+_CHA") |>   # extract only the part like "01_CHA"
  str_replace("_", "")          # remove underscore → "01CHA"

data_SYNGAP$Participant_ID <- data_SYNGAP$Participant_ID |>
  str_extract("^\\d+_CHA") |>   # extract only the part like "01_CHA"
  str_replace("_", "")          # remove underscore → "01CHA"

##renaming age to avoid combination issues
data_hrv <- data_hrv %>%
  rename(Age_HRV = Age)
data_pupil <- data_pupil %>%
  rename(Age_pupil = Age)
data_hrv <- data_hrv %>%
  rename(Participant_ID = Participant.ID)
data_pupil <- data_pupil %>%
  rename(Participant_ID = Participant.ID)

#combine all datasets
data_combined <- data_SYNGAP %>%
  left_join(data_hrv %>% select(-Gender, -Group),
            by = "Participant_ID") %>%
  left_join(data_pupil %>% select(-Gender, -Group),
            by = "Participant_ID")
head(data_combined)


data_combined$Best.Eye <- as.factor(data_combined$Best.Eye)

cols_to_numeric <- c("baseline.NN50", "calming.NN50", "arousing.NN50", 
                     "Lux", "Best.Eye.Valid.Samples.Pre", "Best.Eye.Valid.Samples.Post",
                     "Total.Blinks", "Baseline.Blinks", "Arousing.Blinks", "Calming.Blinks", 
                     "Sample.Rate..Hz.", "Derivative.Threshold..mm.s.")
data_combined[cols_to_numeric] <- lapply(data_combined[cols_to_numeric], function(x) as.numeric(as.character(x)))
str(data_combined) 


physiology_vars <- c(
  "arousing.ibi", "calming.ibi", 
  "arousing.bpm", "calming.bpm",
  "arousing.rmssd", "calming.rmssd", 
  "arousing.pnn50", "calming.pnn50", 
  "arousing.vlf", "calming.vlf", 
  "arousing.lf", "calming.lf",
  "arousing.SD1", "calming.SD1", 
  "arousing.SD2", "calming.SD2")

##filter data combined for SYNGAP only
data_combined_SYNGAP <- data_combined %>%
  filter(Group == "SYNGAP")

##START ANS / PERSISTENCE STATS

#SIB
SIB_ANS_results <- lapply(physiology_vars, function(v) {
  formula <- as.formula(paste(v, "~ SIB_persistence"))
  test <- kruskal.test(formula, data = data_combined_SYNGAP)
  list(variable = v, test = test)
})

SIB_ANS_results

##post hoc test
#BPM
pairwise.wilcox.test(data_combined_SYNGAP$arousing.bpm,
                     data_combined_SYNGAP$SIB_persistence,
                     p.adjust.method = "holm")
#IBI
pairwise.wilcox.test(data_combined_SYNGAP$arousing.ibi,
                     data_combined_SYNGAP$SIB_persistence,
                     p.adjust.method = "holm")

#effect size
library(rstatix)
kruskal_effsize(arousing.bpm ~ SIB_persistence, data = data_combined_SYNGAP)
kruskal_effsize(arousing.ibi ~ SIB_persistence, data = data_combined_SYNGAP)

#direction
#bpm
data_combined_SYNGAP %>%
  group_by(SIB_persistence) %>%
  summarise(
    n = n(),
    median_bpm = median(arousing.bpm, na.rm = TRUE),
    IQR_bpm = IQR(arousing.bpm, na.rm = TRUE)
  )

#ibi
data_combined_SYNGAP %>%
  group_by(SIB_persistence) %>%
  summarise(
    median_ibi = median(arousing.ibi, na.rm = TRUE),
    IQR_ibi = IQR(arousing.ibi, na.rm = TRUE)
  )


#PAG
PAG_ANS_results <- lapply(physiology_vars, function(v) {
  formula <- as.formula(paste(v, "~ PAG_persistence"))
  test <- kruskal.test(formula, data = data_combined_SYNGAP)
  list(variable = v, test = test)
})

PAG_ANS_results


##post hoc tests
#BPM
pairwise.wilcox.test(data_combined_SYNGAP$baseline.pnn50,
                     data_combined_SYNGAP$PAG_persistence,
                     p.adjust.method = "holm")

#effect size
kruskal_effsize(baseline.pnn50 ~ PAG_persistence, data = data_combined_SYNGAP)


#pnn50
data_combined_SYNGAP %>%
  group_by(PAG_persistence) %>%
  summarise(
    n = n(),
    median_bpm = median(baseline.pnn50, na.rm = TRUE),
    IQR_bpm = IQR(baseline.pnn50, na.rm = TRUE)
  )



#STB
STB_ANS_results <- lapply(physiology_vars, function(v) {
  formula <- as.formula(paste(v, "~ STB_persistence"))
  test <- kruskal.test(formula, data = data_combined_SYNGAP)
  list(variable = v, test = test)
})

STB_ANS_results


##post hoc test
#BPM
pairwise.wilcox.test(data_combined_SYNGAP$calming.bpm,
                     data_combined_SYNGAP$STB_persistence,
                     p.adjust.method = "holm")
#IBI
pairwise.wilcox.test(data_combined_SYNGAP$calming.ibi,
                     data_combined_SYNGAP$STB_persistence,
                     p.adjust.method = "holm")

#effect size
kruskal_effsize(calming.bpm ~ STB_persistence, data = data_combined_SYNGAP)
kruskal_effsize(calming.ibi ~ STB_persistence, data = data_combined_SYNGAP)


#bpm
data_combined_SYNGAP %>%
  group_by(STB_persistence) %>%
  summarise(
    n = n(),
    median_bpm = median(calming.bpm, na.rm = TRUE),
    IQR_bpm = IQR(calming.bpm, na.rm = TRUE)
  )

#ibi
data_combined_SYNGAP %>%
  group_by(STB_persistence) %>%
  summarise(
    median_ibi = median(calming.ibi, na.rm = TRUE),
    IQR_ibi = IQR(calming.ibi, na.rm = TRUE)
  )





###Re-do analysis using log-transformed autonomic metrics
vars_to_log <- c(
  "arousing.ibi", "calming.ibi", 
  "arousing.rmssd", "calming.rmssd",
  "arousing.pnn50", "calming.pnn50",
  "arousing.hf", "calming.hf",
  "arousing.sd1", "calming.sd1"
)

#log variables listed above
data_combined <- data_combined %>%
  mutate(
    across(
      all_of(vars_to_log),
      ~ log(.x),
      .names = "log_{.col}"
    )
  )


#select log variables to correlate
log_vars_to_correlate <- c(
  "log_arousing.ibi", "log_calming.ibi",
  "log_arousing.rmssd", "log_calming.rmssd",
  "log_arousing.pnn50", "log_calming.pnn50",
  "log_arousing.hf", "log_calming.hf",
  "log_arousing.sd1", "log_calming.sd1"
)


#SIB - ANOVA
SIB_log_ANOVA_results <- lapply(log_vars_to_correlate, function(v) {
  
  formula <- as.formula(paste(v, "~ SIB_persistence"))
  
  model <- aov(formula, data = data_combined_SYNGAP)
  
  list(
    variable = v,
    anova = summary(model)
  )
})

SIB_log_ANOVA_results



#PAG - ANOVA
PAG_log_ANOVA_results <- lapply(log_vars_to_correlate, function(v) {
  
  formula <- as.formula(paste(v, "~ PAG_persistence"))
  
  model <- aov(formula, data = data_combined_SYNGAP)
  
  list(
    variable = v,
    anova = summary(model)
  )
})

PAG_log_ANOVA_results





#STB - ANOVA
STB_log_ANOVA_results <- lapply(log_vars_to_correlate, function(v) {
  
  formula <- as.formula(paste(v, "~ STB_persistence"))
  
  model <- aov(formula, data = data_combined_SYNGAP)
  
  list(
    variable = v,
    anova = summary(model)
  )
})

STB_log_ANOVA_results
