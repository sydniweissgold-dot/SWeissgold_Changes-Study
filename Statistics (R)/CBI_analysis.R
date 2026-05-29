data <- read.csv("Analysis_1.csv", header=TRUE)

library(tidyverse)

#remove column labelled 'x' (not needed on work computer)
#data = subset(data, select = -c(X))

#remove column 18
data <- data[-18, ]


#order of stats:
##1 - find descriptives - average number of CBI reported by condition
    ##average severity score total 
##2 - types of CBI reported in SYNGAP - highest to lowest 
##3 - average severity score (w / range) for each CBI type

#clean data 
data$Gender <- as.factor(data$Gender)
data$Case.Control <- as.factor(data$Case.Control)

#setting present variables to factors 
CB_present = c('SIB_Present','PAG_Present','VAG_Present','DST_Present','AP_Present','STB_Present','IV_Present',
                'IRC_Present','PIC_Present','ISB_Present','SMR_Present','STL_Present','SIV_Present')
data[,CB_present] <- lapply(data[,CB_present] , factor)

#check that all are factors that need to be 
str(data)


##ANALYSIS###


##demographics##

##sample size
N_total <- nrow(data)
print(N_total)

N_SYNGAP <- sum(data$Case.Control == "SYNGAP")
N_control <- sum(data$Case.Control == "TDC")

##cat prints results, way to print more than 1 at once 
cat("SYNGAP Sample Size:", N_SYNGAP, "Control Sample Size:", N_control)


##gender -1 = male, 2 = female

#total sample genders (SYNGAP & controls)
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
    ##can rank them highest to lowest

#list columns (present) for SYNGAP vs TDC 

##make data frame -> df_SYNGAP has all columns but only SYNGAP participants
##make 2nd data frame -> present_columns -> data frame with only columns that end with present, and only SYNGAP participants 
df_SYNGAP <- subset(data, Case.Control == 'SYNGAP')
present_columns <- names(df_SYNGAP)[endsWith(names(df_SYNGAP), "Present")]
print(present_columns)

#step 1 - create df SYNGAP/present data frame
df_SYNGAP_present <- df_SYNGAP[, present_columns]


#step 2 - go through column by column, if sum > 0, then add column name to a list 
#first change df_SYNGAP_present to integer
df_SYNGAP_present <- sapply(df_SYNGAP_present, as.numeric)

#this is a matrix now, turn it into a dataframe
df_SYNGAP_present_2 <- as.data.frame(df_SYNGAP_present)


##change value from 1/2s to 0/1s 
# Define a function to replace 1s with 0s and 2s with 1s
replace_values <- function(x) {
  ifelse(x == 1, 0, ifelse(x == 2, 1, x))
}
# Apply the function to all columns using lapply
df_SYNGAP_present_2 <- data.frame(lapply(df_SYNGAP_present_2, replace_values))


#start - make an empty list
selected_columns <- c()

# Loop through the columns
for(col_name in names(df_SYNGAP_present_2)) {
  # Check if the sum of the column is greater than the sum of the number of columns (because 
    #as.numeric changed it to 1s and 2s)
  if(sum(df_SYNGAP_present_2[[col_name]]) > 1) {
    # If the sum is greater than 0, add the column name to the list
    selected_columns <- c(selected_columns, col_name)
  }
}
#print list 
selected_columns


##rank them highest to lowest

#from Chat GPT -> calculate sum of each column
column_sums <- sapply(selected_columns, function(col_name) {
  sum(df_SYNGAP_present_2[[col_name]])
  })
#sort column names based on their sums
sorted_column_names <- selected_columns[order(column_sums, decreasing = TRUE)]
#now to add values to column sums (still in order)
column_sums_values <- setNames(column_sums[sorted_column_names], sorted_column_names)
print(column_sums_values)



##percentages 

#% of SYNGAP who displayed each CB
percentage <- colMeans(df_SYNGAP_present_2) * 100
print(percentage)


#% of TDC who displayed CB -- NEED TO DO



##3 - average severity score (w / range) for each CBI type (only looking at SYNGAP, controls not included here - see df_SYNGAP)
library(dplyr)

# Filter columns that end with 'Severity' and calculate the sum for each column

column_names <- names(df_SYNGAP)

##filter columns that end with 'severity'
severity_columns <- grep("Severity$", column_names, value = TRUE)

#sum the columns
mean_of_severity_columns <- df_SYNGAP %>%
  select(all_of(severity_columns)) %>%
  summarise_all(mean, na.rm = TRUE)

print(mean_of_severity_columns)

range(data$Total.CBI.Score)
range(data$Number.of.CB.reported)




