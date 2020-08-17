### SET THE ROOT DIRECTORY

### SET THE ROOT DIRECTORY

library(survival)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)

load("R_out_MSM2017/analysis_calibration_before_after_2010_female.RData")

plot.devel.female <- plot.devel
plot.valid.female <- plot.valid

plot.devel.female <- plot.devel.female + ylim(0,0.2) + ylab("5 year risk")
plot.valid.female <- plot.valid.female + ylim(0,0.2) + ylab("5 year risk")

load("R_out_MSM2017/analysis_calibration_before_after_2010_male.RData")

plot.devel.male <- plot.devel
plot.valid.male <- plot.valid

plot.devel.male <- plot.devel.male + ylim(0,0.2) + ylab("5 year risk")
plot.valid.male <- plot.valid.male + ylim(0,0.2) + ylab("5 year risk")

Figure1_calibration_devel_valid <- ggarrange(plot.devel.female, plot.devel.male, plot.valid.female, plot.valid.male, nrow = 2, ncol = 2,
                                             common.legend = TRUE, legend = "bottom")
Figure1_calibration_devel_valid
ggsave("figures/MSM2017/Figure1_calibration_devel_valid.png",Figure1_calibration_devel_valid,dpi=600,
       height = 7, width = 7)