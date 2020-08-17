### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(survival)
library(foreach)
library(doParallel)
library(ggplot2)
library(reshape2)
library(ggpubr)

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
load("R_out_MSM2017/run_MSM_female_split_fracpoly_V2.RData")

## These are the coefficients of the four models
unweighted.cox.model$coefficients
unweighted.cox.model.timeadj$coefficients
weighted.cox.model$coefficients
weighted.cox.model.timeadj$coefficients


### -----------------------------------------------------------------------------------------------------------------------------
### First going to define four functions that will generate risk scores from each of the models, for a given input dataset
### Also create the datasets which can be used to generate the predicted risks for the development cohort
### -----------------------------------------------------------------------------------------------------------------------------
create_risks_MSM<-function(data.in){
  ## Extract important coefficients and design matrix
  B<-weighted.cox.model$coefficients
  C<-weighted.cox.model$var
  S<-model.matrix(Surv(patid,rep(0,nrow(data.in)), type='right') ~ stat_start365 + Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                    SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                    Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
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
  S<-model.matrix(Surv(patid,rep(0,nrow(data.in)), type='right') ~ stat_start365 + Ethnicity_T0 + Townsend_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                    SBPfrac1_T0 + SBPfrac2_T0 + Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 + CKD345_T0 + 
                    Corticosteroid_use_T0 + T1dia_T0 + T2dia_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + Migraine_T0 + RA_T0 + Sev_men_ill_qrisk_T0 + SBP_var_T0 +
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
pred.data <- data_devel[!duplicated(data_devel$patid),]

## Create a second dataset, where I can manually change stat_start365 to be one for every patient (in original dataset it is zero)
pred.data.2 <- data_devel[!duplicated(data_devel$patid),]
pred.data.2$stat_start365 <- 1

summary(pred.data$stat_start365)
summary(pred.data.2$stat_start365)

### ----------------------------------------------------------------------------------------------------------
### Going to calculate five different predicted risks before preceeding with the comparisons of interest
### ----------------------------------------------------------------------------------------------------------

## Going to add these variables to the pred.data dataframe, as then I can use the CVD_time and CVD_cens variables
## to calculate kaplan meier (observed risks) for calibration plots later

### Risks from MSM, assuming no statin treatment during follow up
pred.data$weighted.risks <- 1 - create_risks_MSM(pred.data)
print("done 3")
### Risks from MSM, assuming no statin treatment during follow up, adjusting for calendar time at baseline
pred.data$weighted.risks.timeadj <- 1 - create_risks_MSM_timeadj(pred.data)
print("done 4")

# ### Finally, risks from the MSM, assuming everyone is on statins during follow up
# pred.data$weighted.risks.on.statin <- 1 - create_risks_MSM(pred.data.2)
# print("done 5")

print("risks calculaled")


### Finally, want to assess calibration of MSM only in patients who arent treated in follow up

## extract patid's for patients who have a statin during follow up
stat_patid <- data_devel$patid[data_devel$stat_start == 1]
stat_patid_dedup <- stat_patid[!duplicated(stat_patid)]

## Create cohort of no statin in follow up and statin in follow up (this datae used to generate risk scores)
pred.data.nostat <- pred.data[!(pred.data$patid %in% stat_patid), ]
pred.data.stat <- pred.data[(pred.data$patid %in% stat_patid), ]
print("str pred.data")
str(pred.data.nostat)

## Do the same for data_devel, which can be used to calculate kaplan meier observed risks 
data_devel.nostat <- data_devel[!(data_devel$patid %in% stat_patid), ]
data_devel.stat <- data_devel[(data_devel$patid %in% stat_patid), ]
print("str data_devel")
str(data_devel.nostat)


# ## Compare different risks of patients when off statin or on statin during follow up
# mean(pred.data$weighted.risks)
# mean(pred.data$weighted.risks.on.statin)


### Create pred risk deciles
## Start by generating a variable which defines which percentile of risk the patients fall into, according to each risk
pred.data.nostat$weighted.centile<-as.integer(cut(pred.data.nostat$weighted.risks, 
                                                  breaks=quantile(pred.data.nostat$weighted.risks,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                                  )))

pred.data.nostat$weighted.timeadj.centile<-as.integer(cut(pred.data.nostat$weighted.risks.timeadj, 
                                                          breaks=quantile(pred.data.nostat$weighted.risks.timeadj,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                                          )))



## Now calculate the average risk within these centiles
cl <- makeCluster(11)
registerDoParallel(11)
predrisk.weighted.centile.nostat<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                           .packages=c("dplyr","mice","tidyr","survival"))
                                   %dopar%{temp<-subset(pred.data.nostat,weighted.centile==input)
                                           return(mean(temp$weighted.risks))
                                           
                                   })
stopCluster(cl)


cl <- makeCluster(11)
registerDoParallel(11)
predrisk.weighted.timeadj.centile.nostat<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                                   .packages=c("dplyr","mice","tidyr","survival"))
                                           %dopar%{temp<-subset(pred.data.nostat,weighted.timeadj.centile==input)
                                                   return(mean(temp$weighted.risks.timeadj))
                                                   
                                           })
stopCluster(cl)

predrisk.weighted.nostat <- unlist(predrisk.weighted.centile.nostat)
predrisk.weighted.timeadj.nostat <- unlist(predrisk.weighted.timeadj.centile.nostat)


### Next get the KM estitaes within these risk groups

## Join the centiles from pred.data.nostat to data_devel
pred.data.nostat.centiles <- pred.data.nostat[,c("patid","weighted.centile","weighted.timeadj.centile")]
data_devel.nostat <- inner_join(data_devel.nostat, pred.data.nostat.centiles, by = "patid")

head(pred.data.nostat.centiles)
head(data_devel.nostat,n=50)

## Define number of days to evaluate observed risk at
numdays <- 1826

## Now calculate km within these centiles
cl <- makeCluster(11)
registerDoParallel(11)
km.weighted.nostat<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                             .packages=c("dplyr","mice","tidyr","survival"))
                     %dopar%{temp.dat<-subset(data_devel.nostat,weighted.centile==input)
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
km.weighted.timeadj.nostat<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                                     .packages=c("dplyr","mice","tidyr","survival"))
                             %dopar%{temp.dat<-subset(data_devel.nostat,weighted.timeadj.centile==input)
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


km.weighted.nostat <- unlist(km.weighted.nostat)
km.weighted.timeadj.nostat <- unlist(km.weighted.timeadj.nostat)



### Now to the calibration plots
# MSM
plot.data.calibration.weighted.nostat<-data.frame(xvals=1:10,km.weighted.nostat,predrisk.weighted.nostat)
colnames(plot.data.calibration.weighted.nostat)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.weighted.nostat<-melt(plot.data.calibration.weighted.nostat,id="xvals")

# Now plot
calibration.weighted.nostat <-ggplot(plot.data.calibration.weighted.nostat) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calendar time not adjusted") + xlab("10th of predicted risk") + ylab("5 year risk")


# MSM timeadj
plot.data.calibration.weighted.timeadj.nostat<-data.frame(xvals=1:10,km.weighted.timeadj.nostat,predrisk.weighted.timeadj.nostat)
colnames(plot.data.calibration.weighted.timeadj.nostat)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.weighted.timeadj.nostat<-melt(plot.data.calibration.weighted.timeadj.nostat,id="xvals")

# Now plot
calibration.weighted.timeadj.nostat <-ggplot(plot.data.calibration.weighted.timeadj.nostat) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calendar time adjusted") + xlab("10th of predicted risk") + ylab("5 year risk")

calibration.weighted.nostat
calibration.weighted.timeadj.nostat

ggsave("figures/MSM2017/female.calibration.weighted.development.nostatcohort.png", calibration.weighted.nostat, dpi = 600)
ggsave("figures/MSM2017/female.calibration.weighted.timeadj.development.nostatcohort.png", calibration.weighted.timeadj.nostat, dpi = 600) 

female.calibration.development.MSM <- ggarrange(calibration.weighted.nostat, calibration.weighted.timeadj.nostat, 
                                               nrow = 1, ncol = 2,
                                               common.legend = TRUE, legend = "bottom")
ggsave("figures/MSM2017/female.calibration.development.MSM.png", female.calibration.development.MSM, 
       dpi = 600, height = 3.5, width = 7) 


# unweighted.data <- data.frame(xvals=1:10,km.unweighted.nostat,predrisk.unweighted.nostat)
# weighted.data <- data.frame(xvals=1:10,km.weighted.nostat,predrisk.weighted.nostat)
# 
# 
# lm(km.unweighted.nostat ~ predrisk.unweighted.nostat, data = unweighted.data)
# 
# lm(km.weighted.nostat ~ predrisk.weighted.nostat, data = weighted.data)


rm(list=setdiff(ls(),list("calibration.weighted.nostat","calibration.weighted.timeadj.nostat","km.weighted.nostat",
                          "km.weighted.timeadj.nostat","predrisk.weighted.nostat","predrisk.weighted.timeadj.nostat")))

save.image("R_out_MSM2017/female_analysis_calibration_development_MSM.RData")
print("saved")