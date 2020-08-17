### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

### This is the third part of deriving the analysis dataset
### We have the baseline data imputed
### and locf/nocb dataset data_anal_predictor_imp
### Need to combine the baseline data with the non baseline data from data_anal_predictor_imp, then do a final locf imputation

### Finally, once we have the final dataset with the time varying variables, we need to extract the variables at baseline, and
### create a new set of variables which are unchanged (baseline variables)

library(mice)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)

load("R_out_MSM2017/derive_analysis_dataset_male_p2.RData")

### Now need to recombine the multiply imputed dataset with the main analysis dataset, and then do LOCF to finish it off
## Get the non baseline variables from data_anal_predictor_imp
not_baseline_data <- data_anal_predictor_imp[duplicated(data_anal_outcomes$patid),]
str(not_baseline_data)

## Make sure both datasets are in the right order, after imputation etc
not_baseline_data <- subset(not_baseline_data, 
                            select = c(patid,pracid,index_date,gender,age_tv,Atrialfib_tv,Atypical_antipsy_med_tv,BMI_tv,
                                       Cholesterol_HDL_Ratio_tv,CKD345_tv,Corticosteroid_use_tv,T1dia_tv,
                                       T2dia_tv,Erec_dysfunc_tv,Ethnicity_tv,FamHis_tv,HIV_tv,Hypertension_tv,
                                       Migraine_tv,RA_tv,Sev_men_ill_qrisk_tv,SBP_tv,SBP_var_tv,Smoking_tv,SLE_tv,Townsend_tv))

baseline_data_imp <- subset(baseline_data_imp, 
                            select = c(patid,pracid,index_date,gender,age_tv,Atrialfib_tv,Atypical_antipsy_med_tv,BMI_tv,
                                       Cholesterol_HDL_Ratio_tv,CKD345_tv,Corticosteroid_use_tv,T1dia_tv,
                                       T2dia_tv,Erec_dysfunc_tv,Ethnicity_tv,FamHis_tv,HIV_tv,Hypertension_tv,
                                       Migraine_tv,RA_tv,Sev_men_ill_qrisk_tv,SBP_tv,SBP_var_tv,Smoking_tv,SLE_tv,Townsend_tv))

## Combine the two
data_anal_predictor_imp_2 <- rbind(baseline_data_imp, not_baseline_data)

## Sort by patid index_date
data_anal_predictor_imp_2 <- arrange(data_anal_predictor_imp_2, patid, index_date)

print("locf imp start")
Sys.time()
## Finnaly do LOCF for the remaining missing variables
data_anal_predictor_imp_3 <- data_anal_predictor_imp_2 %>%
  group_by(patid) %>%
  # First do last observation carried forward
  fill(BMI_tv, SBP_tv, SBP_var_tv, Smoking_tv) %>%
  ungroup() %>%
  as.data.frame()
print("locf imp end")
Sys.time()

str(data_anal_predictor_imp_3)

## Reset the name of data_anal_predictor_imp_3 to data_anal_predictor_imp and clear workspace
data_anal_predictor_imp <- data_anal_predictor_imp_3

rm(list=setdiff(ls(),list("data_anal_predictor_imp","data_anal_outcomes")))
save.image("R_out_MSM2017/derive_analysis_dataset_male_p3.RData")


### STEP 5 ###
### DERIVE TIME INVARIANT DATA FROM THE TIME VARYING VARIABLES (TAKE FIRST ROW) ###
### STEP 5 ###


## Now want to create a dataset that contains only the first entry from each of them, as this will make up the baseline vars
## First create a vector of patid's, shifted, so that i'th element represents patid of the (i+1)th element
patid.shift <- c(0,data_anal_predictor_imp$patid[1:(nrow(data_anal_predictor_imp)-1)])

## Only output rows where patid is not equal to the patid of the previous row (i.e. first row for each patient)
data_time_0 <- data_anal_predictor_imp[!(data_anal_predictor_imp$patid == patid.shift),]


## Want to rename all the columns to have the approiate names, such as _T0
colnames(data_time_0) <- c("patid","pracid","index_date","gender","age_T0","Atrialfib_T0","Atypical_antipsy_med_T0","BMI_T0",
                           "Cholesterol_HDL_Ratio_T0","CKD345_T0","Corticosteroid_use_T0","T1dia_T0",
                           "T2dia_T0","Erec_dysfunc_T0","Ethnicity_T0","FamHis_T0","HIV_T0","Hypertension_T0",
                           "Migraine_T0","RA_T0","Sev_men_ill_qrisk_T0","SBP_T0","SBP_var_T0","Smoking_T0","SLE_T0","Townsend_T0")

## Want to remove gender, pracid, and index date
## These variables will already be in the dataset when we merge this with the other datasets
data_time_0 <- subset(data_time_0, select = -c(gender, pracid, index_date))


### STEP 6 ###
### COMBINE DATASETS INTO ONE DATASET, READY FOR ANALYSIS ###
### STEP 6 ###


## I know have three datasets of interest
## The outcome data: data_anal_outcomes
## The imputed time varying predictor variables: data_anal_predictor_imp
## The variables at baseline (if missing, uses the imputed value)

## I can simply cbind the first two, as they have equivalent number of rows. Then I will do a one to many merge with
## the baseline variables, as this just has one row per patient. Note I want to remove patid from data_anal_predictor_imp,
## so that we don't get the variable patid twice in the dataset
data_anal_predictor_imp <- subset(data_anal_predictor_imp, select = -c(patid))

## Now cbind the two datasets with outcomes and time varying predictor variables
data_anal_outcomes_predictor <- cbind(data_anal_outcomes,data_anal_predictor_imp)

## and merge with the time invariant data
data_anal <- inner_join(data_anal_outcomes_predictor, data_time_0, by = "patid")


## Remove excess stuff and save workspace
rm(list=setdiff(ls(),list("data_anal","data_anal_predictor_imp","data_anal_outcomes","data_time_0")))

save.image("R_out_MSM2017/derive_analysis_dataset_male_p3.RData")
print("image saved")