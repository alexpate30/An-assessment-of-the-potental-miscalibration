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
load("R_out_MSM2017/run_MSM_male_split_fracpoly_V2.RData")

## These are the coefficients of the four models
unweighted.cox.model$coefficients
unweighted.cox.model.timeadj$coefficients
weighted.cox.model$coefficients
weighted.cox.model.timeadj$coefficients


### -----------------------------------------------------------------------------------------------------------------------------
### First going to define four functions that will generate risk scores from each of the models, for a given input dataset
### Also create the datasets which can be used to generate the predicted risks for the development cohort
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

### Risks from standard interval censored cox model
pred.data$unweighted.risks <- 1 - create_risks(pred.data)
print("done 1")
### Risks from standard interval censored cox model, adjusting for calendar time at baseline
pred.data$unweighted.risks.timeadj <- 1 - create_risks_timeadj(pred.data)
print("done 2")


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

print("start calibration")

### ------------------------------------------------------------------------------------------------------------------------
### Now going to look at calibration of the standard cox model, before and after timeadj, to see if it has worked
### ------------------------------------------------------------------------------------------------------------------------


### First do calibration for standard cox in the development cohort
### Already have the predicted values
### Just need to create the KM variables


## First need to merge the centiles into the data_devel dataset (need full data valid dataset in order to fit the survival model)
## as opposed to just baseline data
pred.data.centiles <- pred.data[,c("patid","unweighted.centile","unweighted.timeadj.centile")]

data_devel <- inner_join(data_devel, pred.data.centiles, by = "patid")

## Define numdays that we want to analyse at
numdays <- 1826

## Standard cox model
cl <- makeCluster(11)
registerDoParallel(11)
km.unweighted<-(foreach(input=c(1,2,3,4,5,6,7,8,9,10), .combine=list, .multicombine=TRUE, 
                        .packages=c("dplyr","mice","tidyr","survival"))
                %dopar%{temp.dat<-subset(data_devel,unweighted.centile==input)
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
                        %dopar%{temp.dat<-subset(data_devel,unweighted.timeadj.centile==input)
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
  ggtitle("Calendar time not adjusted") + xlab("10th of predicted risk") + ylab("5 year risk")


# Standard cox timeadj
plot.data.calibration.unweighted.timeadj<-data.frame(xvals=1:10,km.unweighted.timeadj,predrisk.unweighted.timeadj)
colnames(plot.data.calibration.unweighted.timeadj)<-c("xvals","Observed","Predicted")

# Now reshapre into long format
plot.data.calibration.unweighted.timeadj<-melt(plot.data.calibration.unweighted.timeadj,id="xvals")

# Now plot
calibration.unweighted.timeadj <-ggplot(plot.data.calibration.unweighted.timeadj) + 
  geom_point(aes(x=xvals,y=value,shape=variable)) + scale_shape_manual(values=c(1,19)) + 
  ggtitle("Calendar time adjusted") + xlab("10th of predicted risk") + ylab("5 year risk")



ggsave("figures/MSM2017/male.calibration.unweighted.development.png", calibration.unweighted, dpi = 600)
ggsave("figures/MSM2017/male.calibration.unweighted.timeadj.development.png", calibration.unweighted.timeadj, dpi = 600)

male.calibration.development.nonMSM <- ggarrange(calibration.unweighted, calibration.unweighted.timeadj, 
                                                nrow = 1, ncol = 2,
                                                common.legend = TRUE, legend = "bottom")
ggsave("figures/MSM2017/male.calibration.development.nonMSM.png", male.calibration.development.nonMSM, 
       dpi = 600, height = 3.5, width = 7) 

rm(list=setdiff(ls(),list("calibration.unweighted","calibration.unweighted.timeadj","km.unweighted","km.unweighted.timeadj",
                          "predrisk.unweighted","predrisk.unweighted.timeadj")))

save.image("R_out_MSM2017/male_analysis_calibration_development_nonMSM.RData")
print("saved")