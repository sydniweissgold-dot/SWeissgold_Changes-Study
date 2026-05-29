data <- read.csv("Analysis_1.csv", header=TRUE)

library(tidyverse)

#remove column 18
data <- data[-18, ]


data$Gender <- as.factor(data$Gender)
data$Case.Control <- as.factor(data$Case.Control)

CB_present = c('SIB_Present','PAG_Present','VAG_Present','DST_Present','AP_Present','STB_Present','IV_Present',
                'IRC_Present','PIC_Present','ISB_Present','SMR_Present','STL_Present','SIV_Present')
data[,CB_present] <- lapply(data[,CB_present] , factor)

str(data)


##ANALYSIS###


##demographics##

##sample size
N_total <- nrow(data)
print(N_total)

N_SYNGAP <- sum(data$Case.Control == "SYNGAP")
N_control <- sum(data$Case.Control == "TDC")

cat("SYNGAP Sample Size:", N_SYNGAP, "Control Sample Size:", N_control)


##gender -1 = male, 2 = female

#total sample genders
N_gender <- table(data$Gender)
N_male <- N_gender[["1"]]
N_female <- N_gender[["2"]]
cat("Male Sample Size:", N_male, "Female Sample Size:", N_female)

#gender by condition
gender_by_condition <- table(data$Case.Control, data$Gender)
print(gender_by_condition)

gender_by_condition <- data %>%
  group_by(Case.Control, Gender) %>%
  summarise(count = n()) %>%
  ungroup()
print(gender_by_condition)



###ages
age_by_condition <- aggregate(data$Age, by = list(data$Case.Control), FUN = mean)
colnames(age_by_condition) <- c("Condition", "Mean Age")
print(age_by_condition)


#Inferential#


####1 - average number of CBI reported by condition
mean(data[data$Case.Control == "SYNGAP", 'Number.of.CB.reported'])
mean(data[data$Case.Control == "TDC", 'Number.of.CB.reported'])

median(data[data$Case.Control == "SYNGAP", 'Number.of.CB.reported'])
median(data[data$Case.Control == "TDC", 'Number.of.CB.reported'])


#average severity score 
mean(data[data$Case.Control == "SYNGAP", 'Total.CBI.Score'])
mean(data[data$Case.Control == "TDC", 'Total.CBI.Score']) 



####2 - types of CBI reported in SYNGAP

df_SYNGAP <- subset(data, Case.Control == 'SYNGAP')
present_columns <- names(df_SYNGAP)[endsWith(names(df_SYNGAP), "Present")]
print(present_columns)

df_SYNGAP_present <- df_SYNGAP[, present_columns]

df_SYNGAP_present <- sapply(df_SYNGAP_present, as.numeric)

df_SYNGAP_present_2 <- as.data.frame(df_SYNGAP_present)


##change value from 1/2s to 0/1s 
replace_values <- function(x) {
  ifelse(x == 1, 0, ifelse(x == 2, 1, x))
}
df_SYNGAP_present_2 <- data.frame(lapply(df_SYNGAP_present_2, replace_values))


selected_columns <- c()

# Loop through the columns
for(col_name in names(df_SYNGAP_present_2)) {
  if(sum(df_SYNGAP_present_2[[col_name]]) > 1) {
    selected_columns <- c(selected_columns, col_name)
  }
}
selected_columns


##rank them highest to lowest

column_sums <- sapply(selected_columns, function(col_name) {
  sum(df_SYNGAP_present_2[[col_name]])
  })
sorted_column_names <- selected_columns[order(column_sums, decreasing = TRUE)]
column_sums_values <- setNames(column_sums[sorted_column_names], sorted_column_names)
print(column_sums_values)



##percentages 

#% of SYNGAP who displayed each CB
percentage <- colMeans(df_SYNGAP_present_2) * 100
print(percentage)


library(dplyr)

column_names <- names(df_SYNGAP)

severity_columns <- grep("Severity$", column_names, value = TRUE)

mean_of_severity_columns <- df_SYNGAP %>%
  select(all_of(severity_columns)) %>%
  summarise_all(mean, na.rm = TRUE)

print(mean_of_severity_columns)

range(data$Total.CBI.Score)
range(data$Number.of.CB.reported)




