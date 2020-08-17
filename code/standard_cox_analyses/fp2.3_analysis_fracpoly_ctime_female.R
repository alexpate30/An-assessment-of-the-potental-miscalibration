### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(survival)
library(mfp)

### --------------------------------------------------------------------------------------------------------------------------------------------
### This program checks for fractional polynomials of calender time
### --------------------------------------------------------------------------------------------------------------------------------------------

## First load the imputed data, split into development and validation cohorts
load("R_out_MSM2017/standard_cox_datasets_fracpoly_female.RData")

## Check for fractional polynomials of age BMI and SBP, retain all other variables in model, in the development dataset
frac_test_ctime <- mfp(Surv(CVD_time,CVD_cens_R, type='right') ~  fp(c.time_T0) + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + SBPfrac1_T0 + SBPfrac2_T0 +
                         Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 +
                         CKD345_T0 + Corticosteroid_use_T0 + T1dia_T0 + 
                         T2dia_T0 + Ethnicity_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + 
                         Migraine_T0 + RA_T0 + SBP_var_T0 + Sev_men_ill_qrisk_T0 + Smoking_T0 + SLE_T0 + Townsend_T0, 
                                  family=cox, data=data_devel, select=0.05,
                                  keep = c("agefrac1_T0","Atrialfib_T0","Atypical_antipsy_med_T0","BMIfrac1_T0","BMIfrac2_T0","Cholesterol_HDL_Ratio_T0",
                                           "CKD345_T0","Corticosteroid_use_T0","T1dia_T0",
                                           "T2dia_T0","Ethnicity_T0","FamHis_T0","HIV_T0","Hypertension_T0",
                                           "Migraine_T0","RA_T0","Sev_men_ill_qrisk_T0","SBPfrac1_T0","SBPfrac2_T0","SBP_var_T0","Smoking_T0","SLE_T0","Townsend_T0"))

save.image("R_out_MSM2017/standard_cox_fracpoly_ctime_female.RData")
print("image saved")