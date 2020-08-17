### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(mice)
library(foreach)
library(doParallel)
library(dplyr)
library(imputeTS)
library(tidyr)

# -----------------------------------------------------------------------
### This program imputes the baseline data with one stochastic imputation
# -----------------------------------------------------------------------

## First read in the dataset which contains all predictor variables
data_T0 <- read.table("data/CSV2017/analysis_dataset_male_T0.csv",header=TRUE, sep = ",")
data_T0 <- arrange(data_T0, patid)

## Check removed successfully
str(data_T0)

## Make sure variables are in the order we want
data_T0_reducvar <- subset(data_T0,select = c(patid,pracid,index_date_T0,gender,NelsonAalen,CVD_time,CVD_cens,age_T0,Atrialfib_T0,Atypical_antipsy_med_T0,BMI_T0,
                                             Cholesterol_HDL_Ratio_T0,CKD345_T0,Corticosteroid_use_T0,T1dia_T0,
                                             T2dia_T0,Erec_dysfunc_T0,Ethnicity_T0,FamHis_T0,HIV_T0,Hypertension_T0,
                                             Migraine_T0,RA_T0,Sev_men_ill_qrisk_T0,SBP_T0,SBP_var_T0,Smoking_T0,SLE_T0,Townsend_T0))

## Sort out the ethnicity variable, which has some categories which are the same with different names

## Change pakistani to pakistan
data_T0_reducvar$Ethnicity_T0[data_T0_reducvar$Ethnicity_T0 == "pakistani"] <- "pakistan"
## Change oth_asian to asianother
data_T0_reducvar$Ethnicity_T0[data_T0_reducvar$Ethnicity_T0 == "oth_asian"] <- "asianother"
## Set missing things to actually be NA
data_T0_reducvar$Ethnicity_T0[data_T0_reducvar$Ethnicity_T0 %in% c("missing","unclas","other")] <- NA

## drop unused levels
data_T0_reducvar$Ethnicity_T0 <- droplevels(data_T0_reducvar$Ethnicity_T0)


## Next get all variables in the correct format
data_T0_reducvar$Atrialfib_T0 <- as.factor(data_T0_reducvar$Atrialfib_T0)
data_T0_reducvar$Atypical_antipsy_med_T0 <- as.factor(data_T0_reducvar$Atypical_antipsy_med_T0)
data_T0_reducvar$CKD345_T0 <- as.factor(data_T0_reducvar$CKD345_T0)
data_T0_reducvar$Corticosteroid_use_T0 <- as.factor(data_T0_reducvar$Corticosteroid_use_T0)
data_T0_reducvar$T1dia_T0 <- as.factor(data_T0_reducvar$T1dia_T0)
data_T0_reducvar$T2dia_T0 <- as.factor(data_T0_reducvar$T2dia_T0)
data_T0_reducvar$Erec_dysfunc_T0 <- as.factor(data_T0_reducvar$Erec_dysfunc_T0)
data_T0_reducvar$FamHis_T0 <- as.factor(data_T0_reducvar$FamHis_T0)
data_T0_reducvar$HIV_T0 <- as.factor(data_T0_reducvar$HIV_T0)
data_T0_reducvar$Hypertension_T0 <- as.factor(data_T0_reducvar$Hypertension_T0)
data_T0_reducvar$Migraine_T0 <- as.factor(data_T0_reducvar$Migraine_T0)
data_T0_reducvar$RA_T0 <- as.factor(data_T0_reducvar$RA_T0)
data_T0_reducvar$Sev_men_ill_qrisk_T0 <- as.factor(data_T0_reducvar$Sev_men_ill_qrisk_T0)
data_T0_reducvar$Smoking_T0 <- as.factor(data_T0_reducvar$Smoking_T0)
data_T0_reducvar$SLE_T0 <- as.factor(data_T0_reducvar$SLE_T0)
data_T0_reducvar$Townsend_T0 <- as.factor(data_T0_reducvar$Townsend_T0)

## Check structure of dataset
str(data_T0)
str(data_T0_reducvar)


## Now there are some missing NelsonAalen estimates. These are from people who have are censored at a time
## where there are no events, and therefore a probability/cumulative hazard of an event of these time points
## cannot be estimated. There are only 1049 of these people.
## For them we will just use the cumulative hazard estimate associated with the next follow up/event time
## This means sorting by CVD_time, and using nocb imputation
data_T0_reducvar <- arrange(data_T0_reducvar,CVD_time,patid)

## Next the imputation needs to be done
data_T0_reducvar <- data_T0_reducvar %>%
  fill(NelsonAalen, .direction = "up") %>%
  as.data.frame()

## Finally sort by patid again to get in the correct order
data_T0_reducvar <- arrange(data_T0_reducvar,patid)

## Finally remove the unwanted CVD_time var
data_T0_reducvar <- subset(data_T0_reducvar,select = -c(CVD_time))


### DO THE MULTIPLE IMPUTATION ###
### DO THE MULTIPLE IMPUTATION ###
### DO THE MULTIPLE IMPUTATION ###

## Create an empty imputation object, so I can edit the predictor variables matrix, and methods matrix
imp0<-mice(data_T0_reducvar,maxit=0,print=FALSE)
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
pred[,"index_date_T0"]<-0

print ("new predictor matrix")
pred

## Now split in half, becuase it takes about a month to do the imputation when the dataset is full size, compared to one week
set.seed(10)
rand.bern<-rbinom(dim(data_T0_reducvar)[1],1,0.5)
data_T0_reducvar.0 <- data_T0_reducvar[rand.bern==0,]
data_T0_reducvar.1 <- data_T0_reducvar[rand.bern==1,]

## Put two datasets in list format so I can parallelise the imputation (run two processes side by side)
data_T0_reducvar_list <- list(data_T0_reducvar.0,data_T0_reducvar.1)

## Remove old datasets not needed anymore
rm(data_T0_reducvar.0,data_T0_reducvar.1,data_T0_reducvar)

## Now impute
print("stoch imp start")
Sys.time()
cl <- makeCluster(3)
registerDoParallel(3)
data_T0_mice <- (foreach(input=data_T0_reducvar_list, .combine=list, .multicombine=TRUE, 
                          .packages=c("dplyr","mice","tidyr"))
                  %dopar%{mice(input,m=1,meth=meth,pred=pred,seed=500,maxit=10,print=FALSE)
                  })
stopCluster(cl)
print("stoch imp end")
Sys.time()


## Create long format datasets of the imputed data
data_T0_imp.0 <- mice::complete(data_T0_mice[[1]],action='long')
data_T0_imp.1 <- mice::complete(data_T0_mice[[2]],action='long')

## Combine into one dataset
str(data_T0_imp.0)
str(data_T0_imp.1)

print("Combine into one dataset")
data_T0_imp <- rbind(data_T0_imp.0,data_T0_imp.1)
str(data_T0_imp)

## Finally sort the dataset and remove everything that isn't needed
data_anal <- arrange(data_T0_imp,patid)
rm(data_T0_imp.0,data_T0_imp.1,data_T0_mice,data_T0_imp)

## The variable CVD_cens = 1 if a patient is censored, (from SAS), whereas for the R models we want it to = 1, 
## so create the event/censoreing variable as required by R. Also assign the CVD_time variable, time until CVD event
## Must get these variables from data_T0, as I removed them from the dataset before imputation
data_anal$CVD_cens_R = 1 - data_T0$CVD_cens
data_anal$CVD_time = data_T0$CVD_time

## Save image
rm(list=setdiff(ls(),list("data_anal")))
save.image("R_out_MSM2017/impute_baseline_data_male.RData")
print("image saved")