### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(survival)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

load("R_out_MSM2017/analysis_calibration_before_after_2010_timeadj_female.RData")

plot.devel.female <- plot.devel
plot.valid.female <- plot.valid

plot.valid.female <- plot.valid.female + ylim(0,0.1) + ylab("5 year risk")

load("R_out_MSM2017/analysis_calibration_before_after_2010_timeadj_male.RData")

plot.devel.male <- plot.devel
plot.valid.male <- plot.valid

plot.valid.male <- plot.valid.male + ylim(0,0.1) + ylab("5 year risk")

Figure3_calibration_valid_timeadj <- ggarrange(plot.valid.female, plot.valid.male, nrow = 1, ncol = 2,
                                             common.legend = TRUE, legend = "bottom")
Figure3_calibration_valid_timeadj

ggsave("figures/MSM2017/Figure3_calibration_valid_timeadj.png",Figure3_calibration_valid_timeadj,dpi=600,
       height = 3.5, width = 7)