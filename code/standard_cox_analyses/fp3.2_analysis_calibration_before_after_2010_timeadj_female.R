### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY
library(survival)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(reshape2)

### --------------------------------------------------------------------------------------------------------------------------------------------
### This program develops a standard cox model using same criteria as QRISK3 in patients with index dates prior to 2010,
### and then check the calibration in patients with index dates after 2010, but including calender time as a predictor variable,
### the fractional polynomials for calender time are calculated in analysis_fracpoly_ctime.R
### --------------------------------------------------------------------------------------------------------------------------------------------


## First load the imputed data, split into development and validation cohorts
load("R_out_MSM2017/standard_cox_datasets_fracpoly_female.RData")

## Set number of days at which we want to calulcate cvd risks
numdays <- 1826

## Fit the model
fit_model <- coxph(Surv(CVD_time,CVD_cens_R, type='right') ~ c.timefrac1_T0 + c.timefrac2_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                     SBPfrac1_T0 + SBPfrac2_T0 +
                     Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 +
                     CKD345_T0 + Corticosteroid_use_T0 + T1dia_T0 + 
                     T2dia_T0 + Ethnicity_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + 
                     Migraine_T0 + RA_T0 + SBP_var_T0 + Sev_men_ill_qrisk_T0 + Smoking_T0 + SLE_T0 + Townsend_T0, data=data_devel)

## Create function to generate risks for a given model, and set of predictor variables
create_risks<-function(model,datain){
  ## Extract important coefficients and design matrix
  B<-model$coefficients
  C<-model$var
  S<-model.matrix(Surv(CVD_time,CVD_cens_R, type='right') ~ c.timefrac1_T0 + c.timefrac2_T0 + agefrac1_T0 + BMIfrac1_T0 + BMIfrac2_T0 + 
                    SBPfrac1_T0 + SBPfrac2_T0 +
                    Atrialfib_T0 + Atypical_antipsy_med_T0 + Cholesterol_HDL_Ratio_T0 +
                    CKD345_T0 + Corticosteroid_use_T0 + T1dia_T0 + 
                    T2dia_T0 + Ethnicity_T0 + FamHis_T0 + HIV_T0 + Hypertension_T0 + 
                    Migraine_T0 + RA_T0 + SBP_var_T0 + Sev_men_ill_qrisk_T0 + Smoking_T0 + SLE_T0 + Townsend_T0, data=datain)
  S<-S[,-1]
  p<-(S%*%B)[,1]
  lp<-data.frame(lin.pred=p)
  basehaz<-basehaz(model,centered=FALSE)
  basehaz10y<-basehaz[basehaz$time==numdays,]$hazard
  head(basehaz)
  surv<-exp(-basehaz10y*exp(lp))
  return(list("surv" = surv, "lp" = p))
}


## Now generate the survival probabilities for development and validation cohort
res.out.devel <- create_risks(fit_model,data_devel)
res.out.valid <- create_risks(fit_model,data_valid)

## Turn then into risks
risk.all.devel <- 1-res.out.devel$surv
risk.all.valid <- 1-res.out.valid$surv

## Assign correct variable name (it is in dataframe format)
colnames(risk.all.devel)[1] <- "risk.av"
colnames(risk.all.valid)[1] <- "risk.av"

## Summarise risks for development and validation cohort
mean(100*risk.all.devel$risk.av)
mean(100*risk.all.valid$risk.av)

## Assign the other vars req for calculating KM, and patid/pracid
risk.all.devel$patid<-data_devel$patid
risk.all.devel$pracid<-data_devel$pracid
risk.all.devel$CVD_time<-data_devel$CVD_time
risk.all.devel$CVD_cens_R<-data_devel$CVD_cens_R

risk.all.valid$patid<-data_valid$patid
risk.all.valid$pracid<-data_valid$pracid
risk.all.valid$CVD_time<-data_valid$CVD_time
risk.all.valid$CVD_cens_R<-data_valid$CVD_cens_R


## Create a variable that groups per 10th of predicted risk, to be used for calibration plots
risk.all.devel$centile<-as.integer(cut(risk.all.devel$risk.av, 
                                       breaks=quantile(risk.all.devel$risk.av,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                       )))

risk.all.valid$centile<-as.integer(cut(risk.all.valid$risk.av, 
                                      breaks=quantile(risk.all.valid$risk.av,c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
                                      )))


### ---------------------------------------------------
### First going to do calibration of development cohort
### First going to do calibration of development cohort
### ---------------------------------------------------

## Start by calculating kaplan meier (observed) risk for each risk group, defined by centile

## Now calculate the kaplan meier survival for each of these subgroups
cl <- makeCluster(11)
registerDoParallel(11)
km<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                     .packages=c("dplyr","mice","tidyr","survival"))
             %dopar%{temp.dat<-subset(risk.all.devel,centile==input)
                     temp<-survfit(Surv(CVD_time,CVD_cens_R, type='right') ~ 1,
                                   data=temp.dat)
                     # temp$time does not always have a 3653 value, so find a value that it does have
                     numbers<-numdays
                     range<-0
                     arbitrary.numbers<-sort(temp.dat$CVD_time)
                     nearest<-findInterval(numbers, arbitrary.numbers - range)
                     return(1-temp$surv[temp$time==arbitrary.numbers[nearest]]) 
             })
stopCluster(cl)

km.devel<-unlist(km)

## Now need to calculate the average predicted risk per group, defined by centile
cl <- makeCluster(11)
registerDoParallel(11)
predrisk<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                           .packages=c("dplyr","mice","tidyr","survival"))
                   %dopar%{temp<-subset(risk.all.devel,centile==input)
                           return(mean(temp$risk.av))
                           
                   })
stopCluster(cl)

predrisk.devel<-unlist(predrisk)

km.devel
predrisk.devel



## Now to plot them against eachother using ggplot

## Start by creating a dataset using the above vectors
plot.data.devel<-data.frame(xvals=1:10,km.devel,predrisk.devel)

## Name the columns
colnames(plot.data.devel)<-c("percentile","observed","predicted")

# Now reshapre into long format
plot.data.devel<-melt(plot.data.devel,id="percentile")

# Now plot
plot.devel <- ggplot(plot.data.devel) + geom_point(aes(x=percentile,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Female - development cohort") + theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("10th of predicted risk") + scale_x_continuous(breaks=1:10) + ylab("Average 5 year risk")


plot.devel


### ---------------------------------------------------
### Next going to do calibration of validation cohort
### Next going to do calibration of validation cohort
### ---------------------------------------------------


## Now calculate the kaplan meier survival for each of these subgroups
cl <- makeCluster(11)
registerDoParallel(11)
km<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
             .packages=c("dplyr","mice","tidyr","survival"))
     %dopar%{temp.dat<-subset(risk.all.valid,centile==input)
             temp<-survfit(Surv(CVD_time,CVD_cens_R, type='right') ~ 1,
                           data=temp.dat)
             # temp$time does not always have a 3653 value, so find a value that it does have
             numbers<-numdays
             range<-0
             arbitrary.numbers<-sort(temp.dat$CVD_time)
             nearest<-findInterval(numbers, arbitrary.numbers - range)
             return(1-temp$surv[temp$time==arbitrary.numbers[nearest]]) 
     })
stopCluster(cl)

km.valid<-unlist(km)

## Now need to calculate the average predicted risk per group, defined by centile
cl <- makeCluster(11)
registerDoParallel(11)
predrisk<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                   .packages=c("dplyr","mice","tidyr","survival"))
           %dopar%{temp<-subset(risk.all.valid,centile==input)
                   return(mean(temp$risk.av))
                   
           })
stopCluster(cl)

predrisk.valid<-unlist(predrisk)

km.valid
predrisk.valid



### ---------------------------------------------------
### Compare predicted risks for development and validation cohorts
### Compare predicted risks for development and validation cohorts
### ---------------------------------------------------

## Start by creating a dataset using the above vectors
plot.data.valid<-data.frame(xvals=1:10,km.valid,predrisk.valid)

## Name the columns
colnames(plot.data.valid)<-c("percentile","observed","predicted")

# Now reshapre into long format
plot.data.valid<-melt(plot.data.valid,id="percentile")

# Now plot
plot.valid <- ggplot(plot.data.valid) + geom_point(aes(x=percentile,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Female - validation cohort") + theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("10th of predicted risk") + scale_x_continuous(breaks=1:10) + ylab("Average 5 year risk")


plot.valid



## Now create a plot with the predicted risks for development and valid cohorts on the same plot
plot.data.pred.devel.valid<-data.frame(xvals=1:10,predrisk.devel,predrisk.valid)

colnames(plot.data.pred.devel.valid)<-c("percentile","development (post 2010)","validation (pre 2010)")

# Now reshapre into long format
plot.data.pred.devel.valid<-melt(plot.data.pred.devel.valid,id="percentile")

# Now plot
plot.pred.devel.valid <- ggplot(plot.data.pred.devel.valid) + geom_point(aes(x=percentile,y=value,shape=variable)) + 
  scale_shape_manual(values=c(1,19)) + ggtitle("Female - 5 year risk") + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("10th of predicted risk") + scale_x_continuous(breaks=1:10) + ylab("Average 5 year risk")

plot.pred.devel.valid


rm(list=setdiff(ls(), list("plot.valid","plot.devel","plot.pred.devel.valid",
                           "predrisk.devel","predrisk.valid","km.devel","km.valid",
                           "plot.data.valid","plot.data.devel","plot.data.pred.devel.valid")))

save.image("R_out_MSM2017/analysis_calibration_before_after_2010_timeadj_female.RData")


## Want to save the plots
ggsave("figures/MSM2017/standard.cox.pred.devel.valid.timeadj.female.png",plot.pred.devel.valid,dpi=600)
ggsave("figures/MSM2017/standard.cox.calibration.valid.timeadj.female.png",plot.valid,dpi=600)
ggsave("figures/MSM2017/standard.cox.calibration.devel.timeadj.female.png",plot.devel,dpi=600)
