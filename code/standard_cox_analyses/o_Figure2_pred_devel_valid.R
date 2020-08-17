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

plot.pred.devel.valid.female <- plot.pred.devel.valid
plot.pred.devel.valid.female <- plot.pred.devel.valid.female + ylim(0,0.2) + ylab("5 year risk")

load("R_out_MSM2017/analysis_calibration_before_after_2010_male.RData")

plot.pred.devel.valid.male <- plot.pred.devel.valid
plot.pred.devel.valid.male <- plot.pred.devel.valid.male + ylim(0,0.2) + ylab("5 year risk")

Figure2_pred_devel_valid <- ggarrange(plot.pred.devel.valid.female, plot.pred.devel.valid.male, 
                                             nrow = 1, ncol = 2,
                                             common.legend = TRUE, legend = "bottom")

ggsave("figures/MSM2017/Figure2_pred_devel_valid.png",Figure2_pred_devel_valid,dpi=600, 
       height = 3.5, width = 7)