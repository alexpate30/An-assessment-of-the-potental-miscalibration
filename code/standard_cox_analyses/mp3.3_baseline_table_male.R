### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(mice)
library(foreach)
library(doParallel)
library(dplyr)

# --------------------------------------------------------------------------------------------
### This program creates baseline table for the different cohorts (development and validation)
# --------------------------------------------------------------------------------------------

## First read in the dataset which contains all predictor variables
data_T0 <- read.table("data/CSV2017/analysis_dataset_male_T0.csv",header=TRUE, sep = ",")

str(data_T0)

## Create the calendar time variable (number of days from 1st Jan 1998)
data_T0$c.time_T0 <- data_T0$index_date_T0 - 13880

## Check removed successfully
str(data_T0)

## Make sure variables are in the order we want
data_T0 <- subset(data_T0,select = c(patid,pracid,index_date_T0,gender,NelsonAalen,age_T0,Atrialfib_T0,Atypical_antipsy_med_T0,BMI_T0,
                                     Cholesterol_HDL_Ratio_T0,CKD345_T0,Corticosteroid_use_T0,T1dia_T0,
                                     T2dia_T0,Erec_dysfunc_T0,Ethnicity_T0,FamHis_T0,HIV_T0,Hypertension_T0,
                                     Migraine_T0,RA_T0,Sev_men_ill_qrisk_T0,SBP_T0,SBP_var_T0,Smoking_T0,SLE_T0,Townsend_T0,CVD_time,CVD_cens))

## Sort out the ethnicity variable, which has some categories which are the same with different names

## Change pakistani to pakistan
data_T0$Ethnicity_T0[data_T0$Ethnicity_T0 == "pakistani"] <- "pakistan"
## Change oth_asian to asianother
data_T0$Ethnicity_T0[data_T0$Ethnicity_T0 == "oth_asian"] <- "asianother"
## Set missing things to actually be NA
data_T0$Ethnicity_T0[data_T0$Ethnicity_T0 %in% c("missing","unclas")] <- NA

## drop unused levels
data_T0$Ethnicity_T0 <- droplevels(data_T0$Ethnicity_T0)


## Next get all variables in the correct format
data_T0$Atrialfib_T0 <- as.factor(data_T0$Atrialfib_T0)
data_T0$Atypical_antipsy_med_T0 <- as.factor(data_T0$Atypical_antipsy_med_T0)
data_T0$CKD345_T0 <- as.factor(data_T0$CKD345_T0)
data_T0$Corticosteroid_use_T0 <- as.factor(data_T0$Corticosteroid_use_T0)
data_T0$T1dia_T0 <- as.factor(data_T0$T1dia_T0)
data_T0$T2dia_T0 <- as.factor(data_T0$T2dia_T0)
data_T0$FamHis_T0 <- as.factor(data_T0$FamHis_T0)
data_T0$HIV_T0 <- as.factor(data_T0$HIV_T0)
data_T0$Hypertension_T0 <- as.factor(data_T0$Hypertension_T0)
data_T0$Migraine_T0 <- as.factor(data_T0$Migraine_T0)
data_T0$RA_T0 <- as.factor(data_T0$RA_T0)
data_T0$Sev_men_ill_qrisk_T0 <- as.factor(data_T0$Sev_men_ill_qrisk_T0)
data_T0$Smoking_T0 <- as.factor(data_T0$Smoking_T0)
data_T0$SLE_T0 <- as.factor(data_T0$SLE_T0)
data_T0$Townsend_T0 <- as.factor(data_T0$Townsend_T0)


## The variable CVD_cens = 1 if a patient is censored, (from SAS), whereas for the R models we want it to = 1, 
## so create the event/censoreing variable as required by R
data_T0$CVD_cens_R = 1 - data_T0$CVD_cens
data_T0$CVD_cens_R <- as.factor(data_T0$CVD_cens_R )


## Rename all the variabloe to remove the _T0, just for easier reading and coding below
colnames(data_T0) <- c("patid","pracid","index_date_T0","gender","NelsonAalen","age","Atrialfib",
                       "Atypical_antipsy_med","BMI",
                       "Cholesterol_HDL_Ratio","CKD345","Corticosteroid_use","T1dia",
                       "T2dia","Erec_dysfunc","Ethnicity","FamHis","HIV","Hypertension",
                       "Migraine","RA","Sev_men_ill_qrisk_T0","SBP","SBP_var","Smoking","SLE","Townsend","CVD_time","CVD_cens","CVD_cens_R")


## Now create the development and test data, with index dates either side of 2010 (note index_date = 0 is 1960, as day 0 in SAS is 1st Jan 1960)
data_devel<-data_T0[data_T0$index_date_T0 < 18263,]
data_valid<-data_T0[data_T0$index_date_T0 >= 18263,]

N.devel.male <- nrow(data_devel)
N.valid.male <- nrow(data_valid)

### ---------------------------------------------------------------------------
### Now data is in correct form, derive the variables of interest for the table
### ---------------------------------------------------------------------------

## First going to get the mean and standard deviation of all the continuous variables
## First going to get the mean and standard deviation of all the continuous variables
data_devel.mean <- c(mean(data_devel$age),mean(as.numeric(data_devel$SBP),na.rm = TRUE),mean(data_devel$SBP_var,na.rm = TRUE),
                     mean(data_devel$BMI,na.rm = TRUE),mean(data_devel$Cholesterol_HDL_Ratio,na.rm = TRUE))
data_devel.sd <- c(sd(data_devel$age),sd(as.numeric(data_devel$SBP),na.rm = TRUE),sd(data_devel$SBP_var,na.rm = TRUE),
                   sd(data_devel$BMI,na.rm = TRUE),sd(data_devel$Cholesterol_HDL_Ratio,na.rm = TRUE))

data_valid.mean <- c(mean(data_valid$age),mean(as.numeric(data_valid$SBP),na.rm = TRUE),mean(data_valid$SBP_var,na.rm = TRUE),
                     mean(data_valid$BMI,na.rm = TRUE),mean(data_valid$Cholesterol_HDL_Ratio,na.rm = TRUE))
data_valid.sd <- c(sd(data_valid$age),sd(as.numeric(data_valid$SBP),na.rm = TRUE),sd(data_valid$SBP_var,na.rm = TRUE),
                   sd(data_valid$BMI,na.rm = TRUE),sd(data_valid$Cholesterol_HDL_Ratio,na.rm = TRUE))



## Next want to summarise the categorical variables
## Next want to summarise the categorical variables
Smoking.prop <- cbind(100*prop.table(table(data_devel$Smoking)),100*prop.table(table(data_valid$Smoking)))

Townsend.prop <- cbind(100*prop.table(table(data_devel$Townsend)),100*prop.table(table(data_valid$Townsend)))

Ethnicity.prop <- cbind(100*prop.table(table(data_devel$Ethnicity)),100*prop.table(table(data_valid$Ethnicity)))

Hypertension.prop <- c(100*prop.table(table(data_devel$Hypertension))[2],100*prop.table(table(data_valid$Hypertension))[2])

Famhis_lstrict.prop <- c(100*prop.table(table(data_devel$Famhis_lstrict))[2],100*prop.table(table(data_valid$Famhis_lstrict))[2])

T1dia.prop <- c(100*prop.table(table(data_devel$T1dia))[2],100*prop.table(table(data_valid$T1dia))[2])

T2dia.prop <- c(100*prop.table(table(data_devel$T2dia))[2],100*prop.table(table(data_valid$T2dia))[2])

gender.prop <- c(100*prop.table(table(data_devel$gender))[2],100*prop.table(table(data_valid$gender))[2])

Atrialfib.prop <- c(100*prop.table(table(data_devel$Atrialfib))[2],100*prop.table(table(data_valid$Atrialfib))[2])


Atypical_antipsy_med.prop <- c(100*prop.table(table(data_devel$Atypical_antipsy_med))[2],100*prop.table(table(data_valid$Atypical_antipsy_med))[2])


Corticosteroid_use.prop<- c(100*prop.table(table(data_devel$Corticosteroid_use))[2],100*prop.table(table(data_valid$Corticosteroid_use))[2])


CKD345.prop<- c(100*prop.table(table(data_devel$CKD345))[2],100*prop.table(table(data_valid$CKD345))[2])


CKD45.prop<- c(100*prop.table(table(data_devel$CKD45))[2],100*prop.table(table(data_valid$CKD45))[2])


Erec_dysfunc.prop<- c(100*prop.table(table(data_devel$Erec_dysfunc))[2],100*prop.table(table(data_valid$Erec_dysfunc))[2])


HIV.prop <- c(100*prop.table(table(data_devel$HIV))[2],100*prop.table(table(data_valid$HIV))[2])


FamHis.prop <- c(100*prop.table(table(data_devel$FamHis))[2],100*prop.table(table(data_valid$FamHis))[2])


Migraine.prop <- c(100*prop.table(table(data_devel$Migraine))[2],100*prop.table(table(data_valid$Migraine))[2])


RA.prop <- c(100*prop.table(table(data_devel$RA))[2],100*prop.table(table(data_valid$RA))[2])


T1dia.prop <- c(100*prop.table(table(data_devel$T1dia))[2],100*prop.table(table(data_valid$T1dia))[2])


SLE.prop <- c(100*prop.table(table(data_devel$SLE))[2],100*prop.table(table(data_valid$SLE))[2])

Sevmenill.prop <- c(100*prop.table(table(data_devel$Sev_men_ill_qrisk_T0))[2],100*prop.table(table(data_valid$Sev_men_ill_qrisk_T0))[2])

### Now to summarise missingness
### Now to summarise missingness
SBP.miss <- c(100*sum(is.na(data_devel$SBP))/nrow(data_devel),
              100*sum(is.na(data_valid$SBP))/nrow(data_valid))

BMI.miss <- c(100*sum(is.na(data_devel$BMI))/nrow(data_devel),
              100*sum(is.na(data_valid$BMI))/nrow(data_valid))

SBP_var.miss <- c(100*sum(is.na(data_devel$SBP_var))/nrow(data_devel),
                  100*sum(is.na(data_valid$SBP_var))/nrow(data_valid))

Chol.miss <- c(100*sum(is.na(data_devel$Cholesterol_HDL_Ratio))/nrow(data_devel),
               100*sum(is.na(data_valid$Cholesterol_HDL_Ratio))/nrow(data_valid))

Ethnicity.miss <- c(100*sum(is.na(data_devel$Ethnicity))/nrow(data_devel),
                    100*sum(is.na(data_valid$Ethnicity))/nrow(data_valid))

Smoking.miss <- c(100*sum(is.na(data_devel$Smoking))/nrow(data_devel),
                  100*sum(is.na(data_valid$Smoking))/nrow(data_valid))

Townsend.miss <- c(100*sum(is.na(data_devel$Townsend))/nrow(data_devel),
                   100*sum(is.na(data_valid$Townsend))/nrow(data_valid))


## Create continuous table
cont.vars <- data.frame(cbind(paste(round(data_devel.mean,2)," (",round(data_devel.sd,2),")",sep=""),
                              paste(round(data_valid.mean,2)," (",round(data_valid.sd,2),")",sep="")))
colnames(cont.vars) <- c("Development","Validation")
rownames(cont.vars) <- c("Age","SBP","SBP_var","BMI","Chol/HDL")



## Create categorical table
cat.vars <- rbind(Atrialfib.prop,
                  Atypical_antipsy_med.prop,
                  Corticosteroid_use.prop,
                  CKD345.prop,
                  T1dia.prop,
                  T2dia.prop,
                  Ethnicity.prop,
                  FamHis.prop,
                  HIV.prop,
                  Hypertension.prop ,
                  Migraine.prop,
                  RA.prop,
                  Smoking.prop ,
                  SLE.prop,
                  Townsend.prop, 
                  Sevmenill.prop)
colnames(cat.vars) <- c("Development","Validation")

# Create missingness dataset
miss.vars <- rbind(SBP.miss,SBP_var.miss,BMI.miss,Chol.miss,Smoking.miss,Ethnicity.miss,Townsend.miss)
colnames(miss.vars) <- c("Development","Validation")
rownames(miss.vars) <- c("SBP","SBP_var","BMI","Chol/HDL","Smoking","Ethnicity","Townsend")

rm(list=setdiff(ls(),list("N.devel.male","N.valid.male","cont.vars","cat.vars")))

save.image("R_out_MSM2017/baseline_table_male.RData")

