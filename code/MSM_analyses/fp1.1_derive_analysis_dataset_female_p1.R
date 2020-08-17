### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY


## This program loads the outcome variables in analysis formats from the 10 different RData files and combines them
## It then loads all the predictor variables datasets, and combines them also into analysis format
## These are then combined into one dataset, and the first locf (last obs carried forward) imputation is carried out, 
## followed by a nocb (next observation carried backwards)
## This program only deals with the dataset with the time varying variables


## Load relevant packages
library(dplyr)
library(imputeTS)
library(tidyr)
library(mice)

### STEP 1 ###
### LOAD OUTCOME DATA AND PUT IN ANALYSIS FORMAT ###
### STEP 1 ###

## Load resuts in list format, and combine into datasets
load("R_out_MSM2017/derive_outcome_variables1_female.RData")
res_comb1 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables2_female.RData")
res_comb2 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables3_female.RData")
res_comb3 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables4_female.RData")
res_comb4 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables5_female.RData")
res_comb5 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables6_female.RData")
res_comb6 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables7_female.RData")
res_comb7 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables8_female.RData")
res_comb8 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables9_female.RData")
res_comb9 <- do.call("rbind",res_comb)
load("R_out_MSM2017/derive_outcome_variables10_female.RData")
res_comb10 <- do.call("rbind",res_comb)

rm(res_temp,res_comb,i,prog.num)

## Combine into one dataset
data_anal_outcomes <- rbind(res_comb1,res_comb2,res_comb3,res_comb4,res_comb5,
                            res_comb6,res_comb7,res_comb8,res_comb9,res_comb10)
rm(res_comb1, res_comb2, res_comb3, res_comb4, res_comb5, res_comb6, res_comb7, res_comb8, res_comb9, res_comb10)


### STEP 2 ###
### LOAD PREDICTOR VARIABLE DATA AND PUT TIME VARYING VARIABLES INTO ANALYSIS FORMAT ###
### (i.e. different rows for each variable, under the same column name) ###
### STEP 2 ###

## Now load in the predictor variables
data_T0 <- read.table("data/CSV2017/analysis_dataset_female_T0.csv",header=TRUE, sep = ",")
data_T1 <- read.table("data/CSV2017/analysis_dataset_female_T1.csv",header=TRUE, sep = ",")
data_T2 <- read.table("data/CSV2017/analysis_dataset_female_T2.csv",header=TRUE, sep = ",")
data_T3 <- read.table("data/CSV2017/analysis_dataset_female_T3.csv",header=TRUE, sep = ",")
data_T4 <- read.table("data/CSV2017/analysis_dataset_female_T4.csv",header=TRUE, sep = ",")
data_T5 <- read.table("data/CSV2017/analysis_dataset_female_T5.csv",header=TRUE, sep = ",")
data_T6 <- read.table("data/CSV2017/analysis_dataset_female_T6.csv",header=TRUE, sep = ",")
data_T7 <- read.table("data/CSV2017/analysis_dataset_female_T7.csv",header=TRUE, sep = ",")
data_T8 <- read.table("data/CSV2017/analysis_dataset_female_T8.csv",header=TRUE, sep = ",")
data_T9 <- read.table("data/CSV2017/analysis_dataset_female_T9.csv",header=TRUE, sep = ",")


## Next derive the dataset with one line per patient, and the time of censoring/cvd/statin treatment
data_outcomes <- data_T0[,c("patid","CVD_time","CVD_cens","dtcens_num","first_statin_num")]


## With data_outcomes dataset derived, we want to remove some of the extra variables T0 has (the outcome variables) so that
## the T0 - T9 datasets can be combined properly
data_T0 <- data_T0[,c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date_T0","age_T0","Atrialfib_T0",
"Atypical_antipsy_med_T0","BMI_T0","Cholesterol_HDL_Ratio_T0","CKD345_T0","Corticosteroid_use_T0","T1dia_T0",
"T2dia_T0","Erec_dysfunc_T0","Ethnicity_T0","FamHis_T0","HIV_T0","Hypertension_T0","Migraine_T0","RA_T0","Sev_men_ill_qrisk_T0","SBP_T0",
"SBP_var_T0","Smoking_T0","SLE_T0","Townsend_T0")]

## For each dataset relabel the variables to the same thing, tv = time varying
colnames(data_T0) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T1) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T2) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T3) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T4) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T5) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T6) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T7) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T8) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")

colnames(data_T9) <-c("patid","gender","pracid","dtvalid","dtcens","first_statin","first_cvd_all","index_date","age_tv","Atrialfib_tv",
                      "Atypical_antipsy_med_tv","BMI_tv","Cholesterol_HDL_Ratio_tv","CKD345_tv","Corticosteroid_use_tv","T1dia_tv",
                      "T2dia_tv","Erec_dysfunc_tv","Ethnicity_tv","FamHis_tv","HIV_tv","Hypertension_tv","Migraine_tv","RA_tv","Sev_men_ill_qrisk_tv","SBP_tv",
                      "SBP_var_tv","Smoking_tv","SLE_tv","Townsend_tv")



### All datasets now have same column names and same number of rows
### Now if the index date is bigger than patient being censored or having a CVD event we want to remove these observations
## First extract the number for the initial index date (different for each person so is a vector)
index_date_T0 <- data_T0$index_date

## Retain only observations where the difference between the two index dates is smaller than the time until CVD event or censored
## This means no index dates that happen after CVD event or censord are retained
data_T1_red <- data_T1[data_T1$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T2_red <- data_T2[data_T2$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T3_red <- data_T3[data_T3$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T4_red <- data_T4[data_T4$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T5_red <- data_T5[data_T5$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T6_red <- data_T6[data_T6$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T7_red <- data_T7[data_T7$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T8_red <- data_T8[data_T8$index_date - index_date_T0 < data_outcomes$CVD_time, ]
data_T9_red <- data_T9[data_T9$index_date - index_date_T0 < data_outcomes$CVD_time, ]

## Combine into one dataset, I'm thinking it should be the same number of observations as data_anal_outcomes
data_anal_predictor <- rbind(data_T0, data_T1_red, data_T2_red, data_T3_red, data_T4_red, data_T5_red,
                             data_T6_red, data_T7_red, data_T8_red, data_T9_red)

## Need to sort the predictor dataset by patid, then index_date
data_anal_predictor <- arrange(data_anal_predictor, patid, index_date)

## Finally remove patid from data_anal_predictor, so that we don't have patid in the resulting dataset twice
data_anal_predictor <- subset(data_anal_predictor, select = -c(dtvalid, dtcens, first_statin, first_cvd_all))

## Remove all the excess stuff
rm(data_T0, data_T1, data_T2, data_T3, data_T4, data_T5, data_T6, data_T7, data_T8, data_T9, data_T1_red, data_T2_red, data_T3_red, data_T4_red, data_T5_red,
   data_T6_red, data_T7_red, data_T8_red, data_T9_red, index_date_T0)



### STEP 3 ###
## IMPUTATION OF DATA USING LAST OBSERVATION CARRIED FORWARD AND NEXT OBSERVATION CARRIED BACKWARDS ##
## ONLY WORKS IF ATEAST ONE VALUE IS NOT MISSING FOR THAT PATIENT ##
### STEP 3 ###

## Start by doing LOCF and NOCB, this will impute all the different time points as long ##
## as there is at least one time point with a non missing value ##
print("start imp")
Sys.time()

## Want to do last observation carried forward first, followed by next observation carried backwards
## for those still with missing data
data_anal_predictor_imp <- data_anal_predictor %>%
  group_by(patid) %>%
  # First do last observation carried forward
  fill(BMI_tv, Cholesterol_HDL_Ratio_tv, SBP_tv, SBP_var_tv, Smoking_tv) %>%
  # Then do next observation carried backwards, for any that have no non-missing observations before it
  fill(BMI_tv, Cholesterol_HDL_Ratio_tv, SBP_tv, SBP_var_tv, Smoking_tv, .direction = "up") %>%
  ungroup() %>%
  as.data.frame()

print("end imp")
Sys.time()

str(data_anal_predictor)
str(data_anal_predictor_imp)



## Also want to sort out the Ethnicity variable, and also set any missing values to be white
## Next I need to make the dataset suitable to apply the multiple imputation algorithm
## White, missing and unclassified must be combined first
data_anal_predictor_imp$Ethnicity_tv[data_anal_predictor_imp$Ethnicity_tv %in% c("missing","unclas")] <- "white"
## Next change pakistani to pakistan
data_anal_predictor_imp$Ethnicity_tv[data_anal_predictor_imp$Ethnicity_tv == "pakistani"] <- "pakistan"
## Next change oth_asian to asianother
data_anal_predictor_imp$Ethnicity_tv[data_anal_predictor_imp$Ethnicity_tv == "oth_asian"] <- "asianother"

## Drop unused levels and make white the reference level
data_anal_predictor_imp$Ethnicity_tv <- droplevels(data_anal_predictor_imp$Ethnicity_tv)#
data_anal_predictor_imp$Ethnicity_tv <- relevel(data_anal_predictor_imp$Ethnicity_tv, ref = "white")
prop.table(table(data_anal_predictor_imp$Ethnicity_tv))



## Also want to get all variables in the correct format
data_anal_predictor_imp$Atrialfib_tv <- as.factor(data_anal_predictor_imp$Atrialfib_tv)
data_anal_predictor_imp$Atypical_antipsy_med_tv <- as.factor(data_anal_predictor_imp$Atypical_antipsy_med_tv)
data_anal_predictor_imp$CKD345_tv <- as.factor(data_anal_predictor_imp$CKD345_tv)
data_anal_predictor_imp$Corticosteroid_use_tv <- as.factor(data_anal_predictor_imp$Corticosteroid_use_tv)
data_anal_predictor_imp$T1dia_tv <- as.factor(data_anal_predictor_imp$T1dia_tv)
data_anal_predictor_imp$T2dia_tv <- as.factor(data_anal_predictor_imp$T2dia_tv)
data_anal_predictor_imp$Erec_dysfunc_tv <- as.factor(data_anal_predictor_imp$Erec_dysfunc_tv)
data_anal_predictor_imp$FamHis_tv <- as.factor(data_anal_predictor_imp$FamHis_tv)
data_anal_predictor_imp$HIV_tv <- as.factor(data_anal_predictor_imp$HIV_tv)
data_anal_predictor_imp$Hypertension_tv <- as.factor(data_anal_predictor_imp$Hypertension_tv)
data_anal_predictor_imp$Migraine_tv <- as.factor(data_anal_predictor_imp$Migraine_tv)
data_anal_predictor_imp$RA_tv <- as.factor(data_anal_predictor_imp$RA_tv)
data_anal_predictor_imp$Sev_men_ill_qrisk_tv <- as.factor(data_anal_predictor_imp$Sev_men_ill_qrisk_tv)
data_anal_predictor_imp$Smoking_tv <- as.factor(data_anal_predictor_imp$Smoking_tv)
data_anal_predictor_imp$SLE_tv <- as.factor(data_anal_predictor_imp$SLE_tv)
data_anal_predictor_imp$Townsend_tv <- as.factor(data_anal_predictor_imp$Townsend_tv)


rm(list=setdiff(ls(),list("data_anal_predictor_imp","data_anal_outcomes")))
save.image("R_out_MSM2017/derive_analysis_dataset_female_p1.RData")
print("image saved")
