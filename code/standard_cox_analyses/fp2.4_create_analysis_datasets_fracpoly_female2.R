### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

### ----------------------------------------------------------------------------------------------------------
### This program creates adds the required fractional polynomial terms for calendar time to the analysis dataset
### ----------------------------------------------------------------------------------------------------------

## First load the imputed baseline data
load("R_out_MSM2017/standard_cox_datasets_fracpoly_female.RData")

## Development dataset first
data_devel$c.timefrac1_T0 <- ((data_devel$c.time_T0+1)/1000)^-0.5
data_devel$c.timefrac2_T0 <- ((data_devel$c.time_T0+1)/1000)^0.5

## Validation dataset second
data_valid$c.timefrac1_T0 <- ((data_valid$c.time_T0+1)/1000)^-0.5
data_valid$c.timefrac2_T0 <- ((data_valid$c.time_T0+1)/1000)^0.5


## Remove anything excess
rm(list=setdiff(ls(),list("data_devel","data_valid")))

save.image("R_out_MSM2017/standard_cox_datasets_fracpoly_female.RData")
