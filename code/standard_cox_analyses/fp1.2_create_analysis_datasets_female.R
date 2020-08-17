### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

### ----------------------------------------------------------------------------------------------------------
### This program creates the development and validation datasets, splitting the cohort before and after 2010, 
### and curtailing the outcome variable where appropiate
### ----------------------------------------------------------------------------------------------------------

## First load the imputed baseline data
load("R_out_MSM2017/impute_baseline_data_female.RData")

## Create the calendar time variable (number of days from 1st Jan 1998)
data_anal$c.time_T0 <- data_anal$index_date_T0 - 13880

## If the time until event is bigger than 1826 set it to censored, as we are only looking at five year risks (rather than 10)
data_anal$CVD_cens_R[data_anal$CVD_time > 1826] <- 0

## Then change the event time to 1826
data_anal$CVD_time[data_anal$CVD_time > 1826] <- 1826


## Now create the development and test data, with index dates either side of 2010 (note index_date = 0 is 1960, as day 0 in SAS is 1st Jan 1960)
data_devel<-data_anal[data_anal$index_date_T0 < 18263,]
data_valid<-data_anal[data_anal$index_date_T0 >= 18263,]

## For patients in the development cohort, if their follow up goes past 2010, want to censor them there as well
## As in reality, you would not have data on these people past the cut off point

## For these people, set CVD_cens_R to be zero
data_devel$CVD_cens_R[(data_devel$index_date_T0 + data_devel$CVD_time > 18263)] <- 0

## Also, set the CVD_time to be the difference between 18263 (2010-01-01) and their index date
data_devel$CVD_time[(data_devel$index_date_T0 + data_devel$CVD_time > 18263)] <-
  rep(18263, sum(data_devel$index_date_T0 + data_devel$CVD_time > 18263)) - data_devel$index_date_T0[(data_devel$index_date_T0 + data_devel$CVD_time > 18263)]

rm(list=setdiff(ls(),list("data_devel","data_valid")))
save.image("R_out_MSM2017/standard_cox_datasets_female.RData")

## Summarise nuber of events and follow up time, to check if it matches the outcome variables in the MSM dataset (section 3.2)
sum(data_devel$CVD_time)
mean(data_devel$CVD_time)
sum(data_devel$CVD_cens_R)

