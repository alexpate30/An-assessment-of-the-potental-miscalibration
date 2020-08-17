### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY


prog.num <- 6

## Start by loading the T0 dataset, which should also contain the event times, statin times and censoring times
## T0
data_T0 <- read.table("data/CSV2017/analysis_dataset_female_T0.csv",header=TRUE, sep = ",")

## For now all I am interested in is the outcome stuff, so going to remove other vars
colnames(data_T0)
data_outcomes <- data_T0[,c("patid","CVD_time","CVD_cens","dtcens_num","first_statin_num")]


## Start by creating an empty dataset to store results
res_comb <- vector("list",min(200000,nrow(data_outcomes)-200000*(prog.num-1)))

for (i in (200000*(prog.num-1)+1):min(200000*prog.num,nrow(data_outcomes))){
## Create an empty dataset
res_temp <- data.frame(patid = rep(data_outcomes[i,"patid"],10), 
                   tstart = floor(seq(0,9*182.625,182.625)), 
                   tend = floor(seq(182.625,10*182.625,182.625)),
                   event = rep(0,10),
                   censored = rep(0,10),
                   stat_start = rep(0,10),
                   stat_start183 = rep(0,10),
                   stat_start365 = rep(0,10))

## How to edit this minidataframe for a given patient
## First remove any rows that start after a patient is censored or had a CVD event
res_temp <- res_temp[res_temp$tstart < data_outcomes[i,"CVD_time"],]
res_temp

## Also need to set the "tend" to be the value of min_CVD_cens, when the person actually has the event,
## on the final row. However if there is no event or censoring before final follow up, we want to leave it 
## as 1826
res_temp$tend[nrow(res_temp)] <-  min(1826,data_outcomes[i,"CVD_time"])
res_temp

## Now if CVD_comb was lower than dtcens and the event happens prior to 1826,
## assign an event = 1 to the final row.
## If  censoring happens prior to 1826, assign censored = 1 to 
## the final row
## If neither happen prior to 1826, we leave both as zero
if ((data_outcomes[i,"CVD_cens"] == 0) & (data_outcomes[i,"CVD_time"] <= 1826)){
  res_temp$event[nrow(res_temp)] <- 1} else if((data_outcomes[i,"CVD_cens"] == 1) & (data_outcomes[i,"CVD_time"] <= 1826)){
    res_temp$censored[nrow(res_temp)] <- 1}
    

## Now going to create three different stat_start variables
## The first is just if statins have been initiated before the start of the interval
res_temp$stat_start[data_outcomes[i, "first_statin_num"] < res_temp$tstart] <- 1

## The next also required the statin to happen half a year in advance of a CVD event
res_temp$stat_start183[((data_outcomes[i, "first_statin_num"] < res_temp$tstart) & (data_outcomes[i,"CVD_cens"] == 1)) |
                         ((data_outcomes[i, "first_statin_num"] < res_temp$tstart) & (data_outcomes[i,"CVD_cens"] == 0) & 
                            (data_outcomes[i, "first_statin_num"] + 183 < data_outcomes[i, "CVD_time"]))] <- 1

## The next also required the statin to happen half a year in advance of a CVD event
res_temp$stat_start365[((data_outcomes[i, "first_statin_num"] < res_temp$tstart) & (data_outcomes[i,"CVD_cens"] == 1)) |
                         ((data_outcomes[i, "first_statin_num"] < res_temp$tstart) & (data_outcomes[i,"CVD_cens"] == 0) & 
                            (data_outcomes[i, "first_statin_num"] + 365 < data_outcomes[i, "CVD_time"]))] <- 1
                    
## Then need to add res_temp to the master dataset
res_comb[[i]] <- res_temp

if (i%%10000 == 0){print(i)
                   print(Sys.time())}
}

## Remove excess stuff
rm(data_T0,data_outcomes)

save.image("R_out_MSM2017/derive_outcome_variables6_female.RData") 
