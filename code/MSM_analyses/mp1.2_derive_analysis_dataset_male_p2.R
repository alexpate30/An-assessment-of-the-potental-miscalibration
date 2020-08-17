### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

### This is the second part of deriving the analysis dataset
### After doing the locf, there are still lots of variables with missing data, as it is missing at every follow up time
### We therefore extract just the baseline data, and do one stochastic imputation on this
### In part 3, we then do another locf to fill in the remaining gaps
### Note that the weights will not be affected by the stochastically imputed values, as these will not change over time (all imputed
### with the same value)

library(mice)
library(foreach)
library(doParallel)
library(dplyr)
library(imputeTS)
library(tidyr)


## This is to do the stochastic bit of the imputation by splitting the dataset in half and imputing them seperately, 
## just incase doing it as one takes forever
load("R_out_MSM2017/derive_analysis_dataset_male_p1.RData")

## Check removed successfully
str(data_anal_predictor_imp)

## Make sure variables are in the order we want
data_anal_predictor_imp <- subset(data_anal_predictor_imp, 
                                  select = c(patid,pracid,index_date,gender,age_tv,Atrialfib_tv,Atypical_antipsy_med_tv,BMI_tv,
                                             Cholesterol_HDL_Ratio_tv,CKD345_tv,Corticosteroid_use_tv,T1dia_tv,
                                             T2dia_tv,Erec_dysfunc_tv,Ethnicity_tv,FamHis_tv,HIV_tv,Hypertension_tv,
                                             Migraine_tv,RA_tv,Sev_men_ill_qrisk_tv,SBP_tv,SBP_var_tv,Smoking_tv,SLE_tv,Townsend_tv))

### Reduce data_anal_predictor_imp to just the first value for each patient
baseline_data <- data_anal_predictor_imp[!duplicated(data_anal_predictor_imp$patid),]
str(baseline_data)


## Next need to add the outcome variable for the imputation
## The outcome variable here is the time until event/censoring, and the censoring indicator
## It is actually recommended to use the NelsonAalen estimate of survival at the time of censoring, rather than 
## the censoring time
## While in the actual analysis we have interval censored outcome, this variable is just the outcome variable treated
## as if we were to do a standard cox model

## Start by loading the T0 dataset
data_T0 <- read.table("data/CSV2017/analysis_dataset_male_T0.csv",header=TRUE, sep = ",")


## Now there are some missing NelsonAalen estimates. These are from people who have are censored at a time
## where there are no events, and therefore a probability/cumulative hazard of an event of these time points
## cannot be estimated. There are only 1049 of these people.
## For them we will just use the cumulative hazard estimate associated with the next follow up/event time
## This means sorting by CVD_time, and using nocb imputation
data_T0 <- arrange(data_T0,CVD_time,patid)

## Next the imputation needs to be done
data_T0 <- data_T0 %>%
  fill(NelsonAalen, .direction = "up") %>%
  as.data.frame()

## Finally sort by patid again to get in the correct order
data_T0 <- arrange(data_T0,patid)


## Now add the NelsonAalen and CVD_cens variables into the baseline data dataset
baseline_data$NelsonAalen <- data_T0$NelsonAalen
baseline_data$CVD_cens <- data_T0$CVD_cens

## remove data_T0
rm(data_T0)


## Now do a dry run of imputation with no interations to extract imputation criteria to edit
imp0<-mice(baseline_data,maxit=0,print=FALSE)
imp0

### First extract the list of imputation methods for each variables
meth<-imp0$meth
meth

### Next extract the predictor matrix to assign which variables predict which
pred<-imp0$pred

print ("initial predictor matrix")
pred

### Assign which variables that shouldn't be used in the prediction
pred[,"patid"]<-0
pred[,"pracid"]<-0
pred[,"gender"]<-0
pred[,"index_date"]<-0


## Now split in half
set.seed(10)
rand.bern<-rbinom(dim(baseline_data)[1],1,0.5)
baseline_data.0 <- baseline_data[rand.bern==0,]
baseline_data.1 <- baseline_data[rand.bern==1,]

baseline.data.list <- list(baseline_data.0,baseline_data.1)
rm(baseline_data.0,baseline_data.1)

## Now impute seperately
print("stoch imp start")
Sys.time()
cl <- makeCluster(3)
registerDoParallel(3)
baseline.mice <- (foreach(input=baseline.data.list, .combine=list, .multicombine=TRUE, 
                .packages=c("dplyr","mice","tidyr"))
        %dopar%{mice(input,m=1,meth=meth,pred=pred,seed=3000,maxit=10,print=FALSE)
        })
stopCluster(cl)
print("stoch imp end")
Sys.time()

baseline_data_imp.0 <- mice::complete(baseline.mice[[1]],action='long')
baseline_data_imp.1 <- mice::complete(baseline.mice[[2]],action='long')

## Combine into one dataset
str(baseline_data_imp.0)
str(baseline_data_imp.1)

print("Combine into one dataset")
baseline_data_imp <- rbind(baseline_data_imp.0,baseline_data_imp.1)
str(baseline_data_imp)

## Finally sort the dataset and remove everything that isn't needed
baseline_data_imp <- arrange(baseline_data_imp,patid)
rm(baseline_data_imp.0,baseline_data_imp.1,baseline.mice)


## Finally, want to remove the CVD_time and CVD_cens variables from the dataset
## The actual outcome variables we want to use are stored in the data_anal_predictor_imp dataset, which contains interval
## censored outcomes
str(baseline_data_imp)
baseline_data_imp <- subset(baseline_data_imp, 
                        select = c(patid,pracid,index_date,gender,age_tv,Atrialfib_tv,Atypical_antipsy_med_tv,BMI_tv,
                                   Cholesterol_HDL_Ratio_tv,CKD345_tv,Corticosteroid_use_tv,T1dia_tv,
                                   T2dia_tv,Erec_dysfunc_tv,Ethnicity_tv,FamHis_tv,HIV_tv,Hypertension_tv,
                                   Migraine_tv,RA_tv,Sev_men_ill_qrisk_tv,SBP_tv,SBP_var_tv,Smoking_tv,SLE_tv,Townsend_tv))

## Now baseline data imp and data_anal_predictor should have the same number of columns, so can be concatenated
print("check imputed data and data anal predictor have same number of columns")
dim(baseline_data_imp)
dim(data_anal_predictor_imp)


save.image("R_out_MSM2017/derive_analysis_dataset_male_p2.RData")
