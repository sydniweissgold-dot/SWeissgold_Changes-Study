library(dplyr)
library(car)
library(tidyr)
library(stringr)
library(purrr) #for making tidy correlation map
library(broom) #for making tidy correlation map

data_behaviour <- read.csv("data_Q_CBI.csv", header=TRUE)
head(data_behaviour)

data_hrv <- read.csv("0000_CHANGES_HRV_OUTPUT_MANUAL.csv", header=TRUE)
head(data_hrv)

data_pupil <- read.csv("0000_CHANGES_PUPIL_OUTPUT.csv", header=TRUE)
head(data_pupil)

#remove underscore in pupil particpant ID (to match hrv)
data_pupil$Participant.ID <- gsub("_", "", data_pupil$Participant.ID)

#remove year and _ in behaviour data (to match hrv)
data_behaviour$Participant.ID <- data_behaviour$Participant.ID |>
  str_extract("^\\d+_CHA") |>   # extract only the part like "01_CHA"
  str_replace("_", "")          # remove underscore → "01CHA"

##make sure age is in same format in all datasets to avoid combination issues
#was only adding those w/ 'year' (not 'year months' to combined data)
data_behaviour <- data_behaviour %>% mutate(Age = as.character(Age))
data_hrv       <- data_hrv %>% mutate(Age = as.character(Age))
data_pupil     <- data_pupil %>% mutate(Age = as.character(Age))

#combine all datasets
#method is to remove age, gender, group columns from all except hrv dataset
data_combined <- data_behaviour %>%
  full_join(data_hrv %>% select(-Age, -Gender, -Group),
            by = "Participant.ID") %>%
  full_join(data_pupil %>% select(-Age, -Gender, -Group),
            by = "Participant.ID")
# View the combined data
head(data_combined)



##continue formatting data
##remove 'space' after male
data_combined$Gender <- trimws(data_combined$Gender)
data_combined$Gender <- factor(data_combined$Gender, levels = c("Female", "Male"))
levels(data_combined$Gender)

##same for group
data_combined$Group <- trimws(data_combined$Group)
data_combined$Group <- factor(data_combined$Group, levels = c("SYNGAP", "TDC"))
levels(data_combined$Group)


#change factors & numeric
data_combined$Gender <- as.factor(data_combined$Gender)
data_combined$Group <- as.factor(data_combined$Group)
data_combined$ASD <- as.factor(data_combined$ASD)
data_combined$Epilepsy <- as.factor(data_combined$Epilepsy)
data_combined$ADHD <- as.factor(data_combined$ADHD)
data_combined$Best.Eye <- as.factor(data_combined$Best.Eye)

cols_to_numeric <- c("Age", "baseline.NN50", "calming.NN50", "arousing.NN50", 
                     "Lux", "Best.Eye.Valid.Samples.Pre", "Best.Eye.Valid.Samples.Post",
                     "Total.Blinks", "Baseline.Blinks", "Arousing.Blinks", "Calming.Blinks", 
                     "Sample.Rate..Hz.", "Derivative.Threshold..mm.s.")
data_combined[cols_to_numeric] <- lapply(data_combined[cols_to_numeric], function(x) as.numeric(as.character(x)))
str(data_combined) 
str(data_combined$ADHD) #see individual variable

##ASD, ADHD & epilepsy are factors w/ 3 levels:
#the empty columns are level 1, no is level 2, and yes is level 3



###INFERENTIAL STATS###
###REGRESSIONS#### 

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
vars_to_correlate <- c(
  "log_arousing.ibi", "log_calming.ibi",
  "log_arousing.rmssd", "log_calming.rmssd",
  "log_arousing.pnn50", "log_calming.pnn50",
  "log_arousing.hf", "log_calming.hf",
  "log_arousing.sd1", "log_calming.sd1"
)


###create column that combines IBCL & CBCL 
data_combined <- data_combined %>%
  mutate(CBCL.ADHD.Total = coalesce(CBCL.AP, IBCL.ATT))

##subset SYNGAP group
syngap_data <- data_combined %>%
  filter(Group == "SYNGAP")

#example regression:
summary(lm(SrsTotal ~ calming.rmssd, data = syngap_data))



###run regression for SRS 
regression_Srs <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("SrsTotal ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_Srs)

##multiple comparisons
regression_Srs_corr <- regression_Srs %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr")
  )

View(regression_Srs_corr)


##regression for Sensory Profile
regression_SEN <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("SEN.TOT ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_SEN)


#multiple comparisons
regression_SEN_corr <- regression_SEN %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr")
  )

View(regression_SEN_corr)


##regression for CBCL ADHD

regression_ADHD <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("CBCL.ADHD.Total ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_ADHD)

#multiple comparisons
regression_ADHD_corr <- regression_ADHD %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr")
  )

View(regression_ADHD_corr)



##regression for Sleep
regression_CSH <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("CSH.TSD33 ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CSH)

#multiple comparisons
regression_CSH_corr <- regression_CSH %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "fdr")
  )

View(regression_CSH_corr)


##regression for CBI severity

##first - add mean_severity column (code taken from SSBP Behaviour Analysis)
####CREATE CBI SEVERITY AVERAGE !!! 

# Filter columns that end with 'Severity' and calculate the sum for each column
column_names <- names(syngap_data)
##filter columns that end with 'severity'
severity_columns <- grep("Severity$", column_names, value = TRUE)
#average of the columns - i.e. mean IRC_Severity for each participant
mean_of_severity_columns <- syngap_data %>%
  select(all_of(severity_columns)) %>%
  summarise(across(everything(), ~ mean(.x[.x > 0], na.rm = TRUE))) #change so when calculating mean, diving by number of people who are have at least a score of 1 (rather than only NAs)
print(mean_of_severity_columns)

###specifically make Total CBI (w/ average not sum)
participant_means <- syngap_data %>%
  rowwise() %>%
  mutate(mean_severity = mean(c_across(all_of(severity_columns)), na.rm = TRUE)) %>%
  ungroup()

overall_CBI_severity_mean <- mean(participant_means$mean_severity, na.rm = TRUE)

print(overall_CBI_severity_mean)


##create mean severity column for all participants (both groups)
##adapted so mean is calculated: sum of severity scores per participant divided by number of behaviours which have a severity score of 1 or greater
syngap_data <- syngap_data %>%
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



##now run regression
regression_CBI_severity <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("mean_severity ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CBI_severity)


##regression for CBI prevalence
regression_CBI_prevalence <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("Number.of.CB.reported ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CBI_prevalence)

##regression for CBI_SIB
SIB_subset <- syngap_data %>% #first subset for SIB presence
  filter(SIB_Present == 1)

regression_SIB <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("SIB_Severity ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = SIB_subset)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_SIB)

##regression for CBI_PAG
PAG_subset <- syngap_data %>% #first subset for PAG presence
  filter(PAG_Present == 1)

regression_PAG <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("PAG_Severity ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = PAG_subset)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_PAG)

##regression for CBI_STB -- need to subset to only include those who scored on it
STB_subset <- syngap_data %>% #first subset for SIB presence
  filter(STB_Present == 1)

regression_STB <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("STB_Severity ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = STB_subset)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_STB)


##create CBCL's columns -- CHANGE BELOW
syngap_data <- syngap_data %>%
  mutate(
    CBCL.INT.Total = coalesce(CBCL.INT, IBCL.INT),
    CBCL.EXT.Total = coalesce(CBCL.EXT, IBCL.EXT),
    CBCL.ANX.Total = coalesce(CBCL.A.D, IBCL.A.D)
  )

##regression for CBCL_INT
regression_CBCL_INT <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("CBCL.INT.Total ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CBCL_INT)

##regression for CBCL_EXT
regression_CBCL_EXT <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("CBCL.EXT.Total ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CBCL_EXT)

##regression for CBCL anxiety
regression_CBCL_ANX <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("CBCL.ANX.Total ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    # Extract residual degrees of freedom
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CBCL_ANX)

##regression for CBCL depression -- ONLY CBCL (didn't add IBCL)
regression_CBCL_DEP <- map_df(
  vars_to_correlate,
  \(var) {
    formula <- as.formula(
      paste("CBCL.W.D ~", var)   # Hard-coded outcome
    )
    model <- lm(formula, data = syngap_data)
    df_resid <- df.residual(model)
    broom::tidy(model) %>% mutate(predictor_tested = var, df = df_resid)
  }
)

View(regression_CBCL_DEP)




###look at change in pnn50 between calming and arousing w/ behaviour metrics

#calculate change in pNN50
data_combined <- data_combined %>%
  mutate(pnn50_change = calming.pnn50 - arousing.pnn50)

summary(data_combined$pnn50_change)

#check pnn50 change distribution - NON-NORMALLY DISTRIBUTED
hist(data_combined$pnn50_change)
shapiro.test(data_combined$pnn50_change)

#t-test to see if pNN50 changes by group
t.test(pnn50_change ~ Group, data = data_combined)

###need to re-make CBCL INT/EXT, ANX, mean_severity in data_combined(only in data_syngap right now)
data_combined <- data_combined %>%
  mutate(
    CBCL.INT.Total = coalesce(CBCL.INT, IBCL.INT),
    CBCL.EXT.Total = coalesce(CBCL.EXT, IBCL.EXT),
    CBCL.ANX.Total = coalesce(CBCL.A.D, IBCL.A.D)
  )


data_combined <- data_combined %>%
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




#IS PNN50_CHANGE ASSOCIATED WITH BEHAVIOUR METRICS
behaviour_vars <- c(
  "SrsTotal", "CBCL.ADHD.Total", "CBCL.ANX.Total", "CBCL.W.D",
  "CBCL.EXT", "CBCL.INT.Total", "SEN.TOT", "CSH.TSD33",
  "mean_severity", "Number.of.CB.reported"
)

cor_results <- map_df(behaviour_vars, function(var) {
  test <- cor.test(data_combined$pnn50_change,
                   data_combined[[var]],
                   use = "complete.obs",
                   method = "spearman")
  
  tibble(
    variable = var,
    correlation = test$estimate,
    p_value = test$p.value
  )
})

cor_results




