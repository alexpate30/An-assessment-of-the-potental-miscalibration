### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

########################################################################
########################################################################

## In this program I use fractional polynomials for age, BMI and SBP, and
## for the time varying variables where appropiate in calculations of the
## weights

########################################################################
########################################################################


## Load relevant packages
library(dplyr)
library(imputeTS)
library(tidyr)
library(mice)
library(ipw)
library(survival)

### STEP 1 ###
### LOAD DATASET TO GET IN FORM FOR ANALYSIS ###

## Load imputed dataset
load("R_out_MSM2017/derive_analysis_dataset_female_p3.RData")

rm(data_anal_outcomes,data_anal_predictor_imp,data_time_0)

## Data is stored in data_anal


###########
### Start by deriving the calender time variable
###########


## Start by creating the variable which is the calendar time for the index date (starting from 0)
data_anal$c.time_tv <- data_anal$index_date - min(data_anal$index_date)

## Also need to create the c.time_T0 var
## Now I want to extract just the first value and use it as subsequent
## First create a vector of patid's, shifted, so that i'th element represents patid of the (i+1)th element
patid.shift <- c(0,data_anal$patid[1:(nrow(data_anal)-1)])

## Only output rows where patid is not equal to the patid of the previous row (i.e. first row for each patient)
data_baseline <- data_anal[!(data_anal$patid == patid.shift),]

## Reduce to just patid and c.time
data_baseline <- subset(data_baseline, select = c(patid,c.time_tv))

## Want to rename all the columns to have the approiate names, such as _T0
colnames(data_baseline) <- c("patid","c.time_T0")

## Merge with the original dataset
data_anal <- inner_join(data_anal, data_baseline, by = "patid")



###########
### Now calculate the fractional polynomials (same as in non MSM analysis)
###########

## First the variable at time 0
data_anal$c.timefrac1_T0 <- ((data_anal$c.time_T0+1)/1000)^-0.5
data_anal$c.timefrac2_T0 <- ((data_anal$c.time_T0+1)/1000)^0.5

data_anal$agefrac1_T0 <- log(data_anal$age_T0/100)
data_anal$BMIfrac1_T0 <- (data_anal$BMI_T0/10)^-2
data_anal$BMIfrac2_T0 <- log(data_anal$BMI_T0/10)*(data_anal$BMI_T0/10)^-2
data_anal$SBPfrac1_T0 <- (data_anal$SBP_T0/100)^2
data_anal$SBPfrac2_T0 <- (data_anal$SBP_T0/100)^3

## And then the time varying part, only for BMI and SBP
## The time varying variable of age and calendar time is known from time T0, and therefore cannot
## contribute to the weights in the MSM
data_anal$BMIfrac1_tv <- (data_anal$BMI_tv/10)^-2
data_anal$BMIfrac2_tv<- log(data_anal$BMI_tv/10)*(data_anal$BMI_tv/10)^-2
data_anal$SBPfrac1_tv <- (data_anal$SBP_tv/100)^2
data_anal$SBPfrac2_tv <- (data_anal$SBP_tv/100)^3

############
## I want seperate out patients whose first index date is after 2010
## I also want to censor any patients from before 2010, that's data runs over
############

## First just seperate out the two
## Start by creating vectors of the patids, and the first index date for each person
patids_all <- data_anal$patid[!duplicated(data_anal$patid)]
index_dates_all <-data_anal$index_date[!duplicated(data_anal$patid)]

## Now get only the patids where the first index_date < 18263
patids_devel <- patids_all[index_dates_all < 18263]

## Now get the data into development and validation
data_devel <- data_anal[data_anal$patid %in% patids_devel,]
data_valid <- data_anal[!(data_anal$patid %in% patids_devel),]


## Next step is to remove any observations where the index date (time point) is > 18263
data_devel <- data_devel[data_devel$index_date < 18263,]


## Next issue is that we need to know in what situations events happen > 18263, even though the index date
## for that time period was prior to 18263
## Prior to 18263 (then still count it), or after 18263 censor it. 
## Also need to reduce follow up times to be censored at 18263

## Well all the 'tend event times' are actually when the event was, so just need to censor any event times
## after 18263, and bring the follow up times forward

## So if (data_devel$index_date + data_devel$tend - data_devel$tstart > 18263) {means tend is after 18263}
## then we need to set it to zero
data_devel$event[data_devel$index_date + data_devel$tend - data_devel$tstart > 18263] <- 0

## Then also reduce the follow up time so that tend is at 18263
data_devel$tend[data_devel$index_date + data_devel$tend - data_devel$tstart > 18263] <- 
  data_devel$tstart[data_devel$index_date + data_devel$tend - data_devel$tstart > 18263] + 
  (18263 - data_devel$index_date[data_devel$index_date + data_devel$tend - data_devel$tstart > 18263])


### Want to compare follow up times with 3.1
### output the last observation for each patient
test <- data_devel %>%
  group_by(patid) %>%
  slice(n()) %>%
  ungroup
## It does match
sum(test$tend)
mean(test$tend)
## So there is prefect agreement in terms of outcomes between the non MSM and MSM model datasets

### STEP 2 ###
### DO THE ANALYSIS ###
Sys.time()

## Create a variable stat_end365 which indicates if a statin had been used at the end of each period
## Note that for the final period, this variable is irrelevant as it corresponds to the next patient
stat_start_365_shift <- c(data_devel$stat_start365[2:(nrow(data_devel))],0)
data_devel$stat_end365 <- stat_start_365_shift


## I then want to fit a model using only the first 9 rows of each patient (we are not interested in whether the 
## variables at start of the last period cause statin initiaion at the end of the period) ##
## The weights for a row is calculating probability of being on treatment at the end of the period ##
## Given the statin variable we are using corresponds to being on statins at the start of the period ##
## The weight calculated will therefore correspond to the next row ##
## The weight on the first row will always be one, as all patients are not on treatment at the start ##

## I will therefore remove the last element for each patient
data_anal_rem1 <- data_devel %>%
  group_by(patid) %>%
  slice(-n()) %>%
  ungroup()
data_anal_rem1 <- data.frame(data_anal_rem1)

str(data_anal_rem1)

print("calculate weights")
Sys.time()
## Calculate weights on this dataset
weights.obj <- ipwtm(exposure = stat_end365, family = "survival", 
                     numerator = ~ Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                       SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                       Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                       Smoking_T0 + SLE_T0, 
                     denominator = ~ Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                       SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                       Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                       Smoking_T0 + SLE_T0 +
                       Atrialfib_tv + BMIfrac1_tv + BMIfrac2_tv + SBPfrac1_tv + SBPfrac2_tv + CKD345_tv + Hypertension_tv + RA_tv + Smoking_tv, 
                     id = patid, tstart = tstart, timevar = tend,
                     type = "first", data = data_anal_rem1)
head(weights.obj$ipw.weights,n=100)

## Extract the weights themselves in a vector
weights <- weights.obj$ipw.weights

## Now want to combine these weights, with a weight of 1 at the start (first row) for each patient
## Want to start by creating a vector to fill with the weights
weights.all <- rep(0,nrow(data_devel))

## For each element which is the first for by group patid, assign it to being a one
weights.all[!duplicated(data_devel$patid)] <- 1

## Now the number that are == 0 should be the same length as weights
print("check empty slots is same number of entries as weights")
sum(weights.all == 0)
str(weights)

print("assign weights")
Sys.time()
## Now for each element that is equal to zero, assign it as element i of weights
weights.all[weights.all == 0] <- weights
summary(weights)

## After this we should have all the weights in a vector called weights.all

## Can then fit the model
## Fit cox model with no calendar time
print("cox model 1")
Sys.time()
unweighted.cox.model <- coxph(Surv(tstart, tend, event) ~ Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                                SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                                Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                                Smoking_T0 + SLE_T0 + cluster(patid), data = data_devel)
## Fit cox model with calendar time
print("cox model 2")
Sys.time()
unweighted.cox.model.timeadj <- coxph(Surv(tstart, tend, event) ~ Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                                        SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                                        Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                                        Smoking_T0 + SLE_T0 + c.timefrac1_T0 + c.timefrac2_T0 + cluster(patid), data = data_devel)
## Fit MSM model with no calendar time
print("cox model 3")
Sys.time()
weighted.cox.model <- coxph(Surv(tstart, tend, event) ~ stat_start365 + Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                              SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                              Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                              Smoking_T0 + SLE_T0 + cluster(patid), data = data_devel, weights = weights.all)
## Fit MSM model with no calendar time
print("cox model 4")
Sys.time()
weighted.cox.model.timeadj <- coxph(Surv(tstart, tend, event) ~ stat_start365 + Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                                      SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                                      Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                                      Smoking_T0 + SLE_T0 + c.timefrac1_T0 + c.timefrac2_T0 + cluster(patid), data = data_devel, weights = weights.all)

rm(list=setdiff(ls(),list("unweighted.cox.model","unweighted.cox.model.timeadj","weighted.cox.model","weighted.cox.model.timeadj",
                          "weights.all","data_devel","data_valid")))
save.image("R_out_MSM2017/run_MSM_female_split_fracpoly_V2.RData")
print("image saved")