### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(survival)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)

### --------------------------------------------------------------------------------------------------------------------------------------------
### This program takes the MSM models and equivalent time dependent cox models and generates the various risks and comparisons
### Model 1: Interval censored cox (effectively same as standard cox, however we have generated data slightly differently). Important to use
###          this model rather than the original cox model for comparison with MSM, as they are fit on same data
### Model 2: Interval censored cox, adjusting for calender time
### Model 3: Marginal structural model, calculates causal effect of statin initiation in follow up, allows risk to be generated assuming no
###          statin treatment during follow up
### Model 4: Marginal structural model, also adjusting for calender time at baseline
### --------------------------------------------------------------------------------------------------------------------------------------------

## Load the models
load("R_out_MSM2017/run_MSM_male_split_fracpoly_V2.RData")

## These are the coefficients of the four models
unweighted.cox.model$coefficients
unweighted.cox.model.timeadj$coefficients
weighted.cox.model$coefficients
weighted.cox.model.timeadj$coefficients


### -----------------------------------------------------------------------------------------------------------------------------
### First going to define four functions that will generate risk scores from each of the models, for a given input dataset
### Also create the datasets which can be used to generate the predicted risks for the validation cohort
### -----------------------------------------------------------------------------------------------------------------------------

## Define the functions that will generate risks for each of the models
create_risks<-function(data.in){
  ## Extract important coefficients and design matrix
  B<-unweighted.cox.model$coefficients
  C<-unweighted.cox.model$var
  S<-model.matrix(Surv(patid,rep(0,nrow(data.in)), type='right') ~ Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + agefrac2_T0 + BMIfrac1_T0 +
                    SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                    Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + Erec_dysfunc_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                    Smoking_T0 + SLE_T0, data=data.in)
  S<-S[,-1]
  p<-(S%*%B)[,1]
  basehaz<-basehaz(unweighted.cox.model,centered=FALSE)
  basehaz10y<-basehaz[basehaz$time==1826,]$hazard
  head(basehaz)
  surv<-exp(-basehaz10y*exp(p))
  return(surv)
}

create_risks_timeadj<-function(data.in){
  ## Extract important coefficients and design matrix
  B<-unweighted.cox.model.timeadj$coefficients
  C<-unweighted.cox.model.timeadj$var
  S<-model.matrix(Surv(patid,rep(0,nrow(data.in)), type='right') ~ Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + agefrac2_T0 + BMIfrac1_T0 +
                    SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                    Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + Erec_dysfunc_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                    Smoking_T0 + SLE_T0 + c.timefrac1_T0 + c.timefrac2_T0, data=data.in)
  S<-S[,-1]
  p<-(S%*%B)[,1]
  basehaz<-basehaz(unweighted.cox.model.timeadj,centered=FALSE)
  basehaz10y<-basehaz[basehaz$time==1826,]$hazard
  head(basehaz)
  surv<-exp(-basehaz10y*exp(p))
  return(surv)
}

create_risks_MSM<-function(data.in){
  ## Extract important coefficients and design matrix
  B<-weighted.cox.model$coefficients
  C<-weighted.cox.model$var
  S<-model.matrix(Surv(patid,rep(0,nrow(data.in)), type='right') ~ stat_start365 + Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + agefrac2_T0 + BMIfrac1_T0 +
                    SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                    Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + Erec_dysfunc_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                    Smoking_T0 + SLE_T0, data=data.in)
  S<-S[,-1]
  p<-(S%*%B)[,1]
  basehaz<-basehaz(weighted.cox.model,centered=FALSE)
  basehaz10y<-basehaz[basehaz$time==1826,]$hazard
  head(basehaz)
  surv<-exp(-basehaz10y*exp(p))
  return(surv)
}

create_risks_MSM_timeadj<-function(data.in){
  ## Extract important coefficients and design matrix
  B<-weighted.cox.model.timeadj$coefficients
  C<-weighted.cox.model.timeadj$var
  S<-model.matrix(Surv(patid,rep(0,nrow(data.in)), type='right') ~ stat_start365 + Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + agefrac2_T0 + BMIfrac1_T0 +
                    SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                    Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + Erec_dysfunc_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
                    Smoking_T0 + SLE_T0 + c.timefrac1_T0 + c.timefrac2_T0, data=data.in)
  S<-S[,-1]
  p<-(S%*%B)[,1]
  basehaz<-basehaz(weighted.cox.model.timeadj,centered=FALSE)
  basehaz10y<-basehaz[basehaz$time==1826,]$hazard
  head(basehaz)
  surv<-exp(-basehaz10y*exp(p))
  return(surv)
}


## Only want the first row of information for each patient to generate risks
pred.data <- data_valid[!duplicated(data_valid$patid),]

## Create a second dataset, where I can manually change stat_start365 to be one for every patient (in original dataset it is zero)
pred.data.2 <- data_valid[!duplicated(data_valid$patid),]
pred.data.2$stat_start365 <- 1

summary(pred.data$stat_start365)
summary(pred.data.2$stat_start365)

### ----------------------------------------------------------------------------------------------------------
### Going to calculate five different predicted risks before preceeding with the comparisons of interest
### ----------------------------------------------------------------------------------------------------------

## Going to add these variables to the pred.data dataframe, as then I can use the CVD_time and CVD_cens variables
## to calculate kaplan meier (observed risks) for calibration plots later

### Risks from standard interval censored cox model
pred.data$unweighted.risks <- 1 - create_risks(pred.data)
print("done 1")
### Risks from standard interval censored cox model, adjusting for calendar time at baseline
pred.data$unweighted.risks.timeadj <- 1 - create_risks_timeadj(pred.data)
print("done 2")
### Risks from MSM, assuming no statin treatment during follow up
pred.data$weighted.risks <- 1 - create_risks_MSM(pred.data)
print("done 3")
### Risks from MSM, assuming no statin treatment during follow up, adjusting for calendar time at baseline
pred.data$weighted.risks.timeadj <- 1 - create_risks_MSM_timeadj(pred.data)
print("done 4")

### Finally, risks from the MSM, assuming everyone is on statins during follow up
pred.data$weighted.risks.on.statin <- 1 - create_risks_MSM(pred.data.2)
print("done 5")

print("risks calculaled")

### ----------------------------------------------------------------------------------------------------------
### Compare when everyone is on statins, vs not on statins, as a QC of the MSM
### ----------------------------------------------------------------------------------------------------------

##  Compare risk off statin vs risk on
100*mean(pred.data$weighted.risks)
100*mean(pred.data$weighted.risks.on.statin)

## Calculate the average relative risk
mean(pred.data$weighted.risks.on.statin)/mean(pred.data$weighted.risks)


### --------------------------------------------------------------------------------------------------------------
### Next quantify the increase in risk we get from the standard cox model, when using the MSM (stat_start365 = 0)
### which assumes no statin use during follow up
### --------------------------------------------------------------------------------------------------------------
mean(pred.data$weighted.risks)
mean(pred.data$unweighted.risks)

ratio.no.statin <- mean(pred.data$weighted.risks)/mean(pred.data$unweighted.risks)
ratio.no.statin

## An increase of %
100*(mean(pred.data$weighted.risks)/mean(pred.data$unweighted.risks) - 1)

### --------------------------------------------------------------------------------------------------------------
### Next quantify the reduction in risk we see after introducing the calender time variable
### --------------------------------------------------------------------------------------------------------------

print("quantify risk reductions")

## Start by comparing risks before and after timeadj and producing a ratio
100*mean(pred.data$unweighted.risks)
100*mean(pred.data$unweighted.risks.timeadj)

## Get proportion of the original risks that the time adjusted risks are
predrisk.comparison.unweighted <- 100*(mean(pred.data$unweighted.risks)-mean(pred.data$unweighted.risks.timeadj))/mean(pred.data$unweighted.risks)
predrisk.comparison.unweighted

## Start by generating a variable which defines which percentile of risk the patients fall into, according to each risk
pred.data$unweighted.centile<-as.integer(cut(pred.data$unweighted.risks, 
                                             breaks=quantile(pred.data$unweighted.risks,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                             )))

pred.data$unweighted.timeadj.centile<-as.integer(cut(pred.data$unweighted.risks.timeadj, 
                                                     breaks=quantile(pred.data$unweighted.risks.timeadj,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                                     )))

## Now calculate the average risk within these centiles
cl <- makeCluster(11)
registerDoParallel(11)
predrisk.unweighted.centile<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                      .packages=c("dplyr","mice","tidyr","survival"))
                              %dopar%{temp<-subset(pred.data,unweighted.centile==input)
                                      return(mean(temp$unweighted.risks))
                                      
                              })
stopCluster(cl)


cl <- makeCluster(11)
registerDoParallel(11)
predrisk.unweighted.timeadj.centile<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                              .packages=c("dplyr","mice","tidyr","survival"))
                                      %dopar%{temp<-subset(pred.data,unweighted.timeadj.centile==input)
                                              return(mean(temp$unweighted.risks.timeadj))
                                              
                                      })
stopCluster(cl)


predrisk.unweighted <- unlist(predrisk.unweighted.centile)
predrisk.unweighted.timeadj <- unlist(predrisk.unweighted.timeadj.centile)

# Now make a data frame with everything I want to plot
plot.data.unweighted<-data.frame(xvals=1:10,predrisk.unweighted,predrisk.unweighted.timeadj)
colnames(plot.data.unweighted)<-c("xvals","Not adjusting for time","Adjusting for time")

# Now reshapre into long format
plot.data.unweighted<-melt(plot.data.unweighted,id="xvals")

# Now plot
predrisk.comparison.unweighted.centile <- ggplot(plot.data.unweighted) +
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Comparison of risks after adjusting for calendar time in standard cox")
predrisk.comparison.unweighted.centile


### ------------------------------------------------------------------------------------------------------------------------
### Next quantify the reduction in risk we see after introducing the calendar time variable, in the MSM Setting (no statins)
### ------------------------------------------------------------------------------------------------------------------------

## Start by comparing risks before and after timeadj and producing a ratio
100*mean(pred.data$weighted.risks)
100*mean(pred.data$weighted.risks.timeadj)

## Get proportion of the original risks that the time adjusted risks are
predrisk.comparison.weighted <- 100*(mean(pred.data$weighted.risks)-mean(pred.data$weighted.risks.timeadj))/mean(pred.data$weighted.risks)
predrisk.comparison.weighted

## Start by generating a variable which defines which percentile of risk the patients fall into, according to each risk
pred.data$weighted.centile<-as.integer(cut(pred.data$weighted.risks, 
                                             breaks=quantile(pred.data$weighted.risks,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                             )))

pred.data$weighted.timeadj.centile<-as.integer(cut(pred.data$weighted.risks.timeadj, 
                                                     breaks=quantile(pred.data$weighted.risks.timeadj,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                                     )))

## Now calculate the average risk within these centiles
cl <- makeCluster(11)
registerDoParallel(11)
predrisk.weighted.centile<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                      .packages=c("dplyr","mice","tidyr","survival"))
                              %dopar%{temp<-subset(pred.data,weighted.centile==input)
                                      return(mean(temp$weighted.risks))
                                      
                              })
stopCluster(cl)


cl <- makeCluster(11)
registerDoParallel(11)
predrisk.weighted.timeadj.centile<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                              .packages=c("dplyr","mice","tidyr","survival"))
                                      %dopar%{temp<-subset(pred.data,weighted.timeadj.centile==input)
                                              return(mean(temp$weighted.risks.timeadj))
                                              
                                      })
stopCluster(cl)


predrisk.weighted <- unlist(predrisk.weighted.centile)
predrisk.weighted.timeadj <- unlist(predrisk.weighted.timeadj.centile)

# Now make a data frame with everything I want to plot
plot.data.weighted<-data.frame(xvals=1:10,predrisk.weighted,predrisk.weighted.timeadj)
colnames(plot.data.weighted)<-c("xvals","Not adjusting for time","Adjusting for time")

# Now reshapre into long format
plot.data.weighted<-melt(plot.data.weighted,id="xvals")

# Now plot
predrisk.comparison.weighted.centile <- ggplot(plot.data.weighted) +
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Comparison of risks after adjusting for calendar time in standard cox")
predrisk.comparison.weighted.centile



##################################################################################################################
### --------------------------------------------------------------------------------------------------------------
### CALIBRATION STUFF
### CALIBRATION STUFF
### CALIBRATION STUFF
### --------------------------------------------------------------------------------------------------------------
##################################################################################################################

print("start calibration")

### ------------------------------------------------------------------------------------------------------------------------
### Now going to look at calibration of the standard cox model, before and after timeadj, to see if it has worked
### ------------------------------------------------------------------------------------------------------------------------


### First do calibration for standard cox in the validation cohort
### Already have the predicted values
### Just need to create the KM variables


## First need to merge the centiles into the data_valid dataset (need full data valid dataset in order to fit the survival model)
## as opposed to just baseline data
pred.data.centiles <- pred.data[,c("patid","unweighted.centile","unweighted.timeadj.centile",
                                   "weighted.centile","weighted.timeadj.centile")]

data_valid <- inner_join(data_valid, pred.data.centiles, by = "patid")

## Define numdays that we want to analyse at
numdays <- 1826

## Standard cox model
cl <- makeCluster(11)
registerDoParallel(11)
km.unweighted<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                        .packages=c("dplyr","mice","tidyr","survival"))
                %dopar%{temp.dat<-subset(data_valid,unweighted.centile==input)
                        temp<-survfit(Surv(tstart,tend, event) ~ 1,
                                      data=temp.dat)
                        # temp$time does not always have a 3653 value, so find a value that it does have
                        numbers<-numdays
                        range<-0
                        arbitrary.numbers<-sort(temp.dat$tend)
                        nearest<-findInterval(numbers, arbitrary.numbers - range)
                        return(1-temp$surv[temp$time==arbitrary.numbers[nearest]]) 
                })
stopCluster(cl)


## Standard cox adjusting for calendar time
cl <- makeCluster(11)
registerDoParallel(11)
km.unweighted.timeadj<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                .packages=c("dplyr","mice","tidyr","survival"))
                        %dopar%{temp.dat<-subset(data_valid,unweighted.timeadj.centile==input)
                                temp<-survfit(Surv(tstart,tend, event) ~ 1,
                                              data=temp.dat)
                                # temp$time does not always have a 3653 value, so find a value that it does have
                                numbers<-numdays
                                range<-0
                                arbitrary.numbers<-sort(temp.dat$tend)
                                nearest<-findInterval(numbers, arbitrary.numbers - range)
                                return(1-temp$surv[temp$time==arbitrary.numbers[nearest]]) 
                        })
stopCluster(cl)


km.unweighted <- unlist(km.unweighted)
km.unweighted.timeadj <- unlist(km.unweighted.timeadj)


### Calibration plots

# Standard cox
plot.data.calibration.unweighted<-data.frame(xvals=1:10,km.unweighted,predrisk.unweighted)
colnames(plot.data.calibration.unweighted)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.unweighted<-melt(plot.data.calibration.unweighted,id="xvals")

# Now plot
calibration.unweighted <-ggplot(plot.data.calibration.unweighted) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calibration (validation cohort) - standard cox")


# Standard cox timeadj
plot.data.calibration.unweighted.timeadj<-data.frame(xvals=1:10,km.unweighted.timeadj,predrisk.unweighted.timeadj)
colnames(plot.data.calibration.unweighted.timeadj)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.unweighted.timeadj<-melt(plot.data.calibration.unweighted.timeadj,id="xvals")

# Now plot
calibration.unweighted.timeadj <-ggplot(plot.data.calibration.unweighted.timeadj) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calibration (validation cohort) - Standard cox, calendar time adjusted")



## MSM
cl <- makeCluster(11)
registerDoParallel(11)
km.weighted<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                      .packages=c("dplyr","mice","tidyr","survival"))
              %dopar%{temp.dat<-subset(data_valid,weighted.centile==input)
                      temp<-survfit(Surv(tstart,tend, event) ~ 1,
                                    data=temp.dat)
                      # temp$time does not always have a 3653 value, so find a value that it does have
                      numbers<-numdays
                      range<-0
                      arbitrary.numbers<-sort(temp.dat$tend)
                      nearest<-findInterval(numbers, arbitrary.numbers - range)
                      return(1-temp$surv[temp$time==arbitrary.numbers[nearest]]) 
              })
stopCluster(cl)


## Standard cox adjusting for calendar time
cl <- makeCluster(11)
registerDoParallel(11)
km.weighted.timeadj<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                              .packages=c("dplyr","mice","tidyr","survival"))
                      %dopar%{temp.dat<-subset(data_valid,weighted.timeadj.centile==input)
                              temp<-survfit(Surv(tstart,tend, event) ~ 1,
                                            data=temp.dat)
                              # temp$time does not always have a 3653 value, so find a value that it does have
                              numbers<-numdays
                              range<-0
                              arbitrary.numbers<-sort(temp.dat$tend)
                              nearest<-findInterval(numbers, arbitrary.numbers - range)
                              return(1-temp$surv[temp$time==arbitrary.numbers[nearest]]) 
                      })
stopCluster(cl)


km.weighted <- unlist(km.weighted)
km.weighted.timeadj <- unlist(km.weighted.timeadj)

# MSM
plot.data.calibration.weighted<-data.frame(xvals=1:10,km.weighted,predrisk.weighted)
colnames(plot.data.calibration.weighted)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.weighted<-melt(plot.data.calibration.weighted,id="xvals")

# Now plot
calibration.weighted <-ggplot(plot.data.calibration.weighted) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calibration (validation cohort) - MSM")


# MSM timeadj
plot.data.calibration.weighted.timeadj<-data.frame(xvals=1:10,km.weighted.timeadj,predrisk.weighted.timeadj)
colnames(plot.data.calibration.weighted.timeadj)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.weighted.timeadj<-melt(plot.data.calibration.weighted.timeadj,id="xvals")

# Now plot
calibration.weighted.timeadj <-ggplot(plot.data.calibration.weighted.timeadj) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calibration (validation cohort) - MSM, calendar time adjusted")


ggsave("figures/MSM2017/male.calibration.unweighted.fracpoly_V2.png", calibration.unweighted, dpi = 600)
ggsave("figures/MSM2017/male.calibration.unweighted.timeadj.fracpoly_V2.png", calibration.unweighted.timeadj, dpi = 600)
ggsave("figures/MSM2017/male.calibration.weighted.fracpoly_V2.png", calibration.weighted, dpi = 600)
ggsave("figures/MSM2017/male.calibration.weighted.timeadj.fracpoly_V2.png", calibration.weighted.timeadj, dpi = 600) 

rm(data_devel, data_valid)
save.image("R_out_MSM2017/analysis_MSM_male_split_fracpoly_V2.RData")
print("saved")