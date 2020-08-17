### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

### ----------------------------------------------------------------------------------------------------------
### This program creates adds the required fractional polynomial terms to the analysis dataset
### ----------------------------------------------------------------------------------------------------------

## First load the imputed baseline data
load("R_out_MSM2017/standard_cox_datasets_male.RData")

## Development dataset first
data_devel$agefrac1_T0 <- (data_devel$age_T0/100)^-1
data_devel$agefrac2_T0 <- (data_devel$age_T0/100)^3
data_devel$BMIfrac1_T0 <- (data_devel$BMI_T0/10)
data_devel$SBPfrac1_T0 <- (data_devel$SBP_T0/100)^3
data_devel$SBPfrac2_T0 <- ((data_devel$SBP_T0/100)^3)*log(data_devel$SBP_T0/100)

## Validation dataset second
data_valid$agefrac1_T0 <- (data_valid$age_T0/100)^-1
data_valid$agefrac2_T0 <- (data_valid$age_T0/100)^3
data_valid$BMIfrac1_T0 <- (data_valid$BMI_T0/10)
data_valid$SBPfrac1_T0 <- (data_valid$SBP_T0/100)^3
data_valid$SBPfrac2_T0 <- ((data_valid$SBP_T0/100)^3)*log(data_valid$SBP_T0/100)

## Remove anything excess
rm(list=setdiff(ls(),list("data_devel","data_valid")))

save.image("R_out_MSM2017/standard_cox_datasets_fracpoly_male.RData")


