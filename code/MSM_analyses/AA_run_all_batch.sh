LOAD R
module load apps/binapps/sas/9.4
module load apps/gcc/R/3.4.2

cd "ROOT DIRECTORY/code"
nohup Rscript fp0.1_derive_outcome_variables1_female.R > fp0.1_derive_outcome_variables1_female.out
nohup Rscript fp0.2_derive_outcome_variables2_female.R > fp0.2_derive_outcome_variables2_female.out
nohup Rscript fp0.3_derive_outcome_variables3_female.R > fp0.3_derive_outcome_variables3_female.out
nohup Rscript fp0.4_derive_outcome_variables4_female.R > fp0.4_derive_outcome_variables4_female.out
nohup Rscript fp0.5_derive_outcome_variables5_female.R > fp0.5_derive_outcome_variables5_female.out
nohup Rscript fp0.6_derive_outcome_variables6_female.R > fp0.6_derive_outcome_variables6_female.out
nohup Rscript fp0.7_derive_outcome_variables7_female.R > fp0.7_derive_outcome_variables7_female.out
nohup Rscript fp0.8_derive_outcome_variables8_female.R > fp0.8_derive_outcome_variables8_female.out
nohup Rscript fp0.9_derive_outcome_variables9_female.R > fp0.9_derive_outcome_variables9_female.out
nohup Rscript fp0.10_derive_outcome_variables10_female.R > fp0.10_derive_outcome_variables10_female.out
nohup Rscript fp1.1_derive_analysis_dataset_female_p1.R > fp1.1_derive_analysis_dataset_female_p1.out
nohup Rscript fp1.2_derive_analysis_dataset_female_p2.R > fp1.2_derive_analysis_dataset_female_p2.out
nohup Rscript fp1.3_derive_analysis_dataset_female_p3.R > fp1.3_derive_analysis_dataset_female_p3.out
nohup Rscript fp2.1_run_MSM_female_split.R > fp2.1_run_MSM_female_split.out 
nohup Rscript fp2.2_analysis_MSM_generate_risks_comparisons_female.R > fp2.2_analysis_MSM_generate_risks_comparisons_female.out
nohup Rscript fp2.3_analysis_calibration_development_MSM.R > fp2.3_analysis_calibration_development_MSM.
nohup Rscript fp2.3_analysis_calibration_development_nonMSM.R > fp2.3_analysis_calibration_development_nonMSM.out
nohup Rscript fp2.3_analysis_calibration_validation_MSM.R > fp2.3_analysis_calibration_validation_MSM.
nohup Rscript fp2.3_analysis_calibration_validation_nonMSM.R > fp2.3_analysis_calibration_validation_nonMSM.out
nohup Rscript mp0.1_derive_outcome_variables1_male.R > mp0.1_derive_outcome_variables1_male.out
nohup Rscript mp0.2_derive_outcome_variables2_male.R > mp0.2_derive_outcome_variables2_male.out
nohup Rscript mp0.3_derive_outcome_variables3_male.R > mp0.3_derive_outcome_variables3_male.out
nohup Rscript mp0.4_derive_outcome_variables4_male.R > mp0.4_derive_outcome_variables4_male.out
nohup Rscript mp0.5_derive_outcome_variables5_male.R > mp0.5_derive_outcome_variables5_male.out
nohup Rscript mp0.6_derive_outcome_variables6_male.R > mp0.6_derive_outcome_variables6_male.out
nohup Rscript mp0.7_derive_outcome_variables7_male.R > mp0.7_derive_outcome_variables7_male.out
nohup Rscript mp0.8_derive_outcome_variables8_male.R > mp0.8_derive_outcome_variables8_male.out
nohup Rscript mp0.9_derive_outcome_variables9_male.R > mp0.9_derive_outcome_variables9_male.out
nohup Rscript mp0.10_derive_outcome_variables10_male.R > mp0.10_derive_outcome_variables10_male.
nohup Rscript mp1.1_derive_analysis_dataset_male_p1.R > mp1.1_derive_analysis_dataset_male_p1.out
nohup Rscript mp1.2_derive_analysis_dataset_male_p2.R > mp1.2_derive_analysis_dataset_male_p2.out
nohup Rscript mp1.3_derive_analysis_dataset_male_p3.R > mp1.3_derive_analysis_dataset_male_p3.out
nohup Rscript mp2.1_run_MSM_male_split.R > mp2.1_run_MSM_male_split.out
nohup Rscript mp2.2_analysis_MSM_generate_risks_comparisons_male.R > mp2.2_analysis_MSM_generate_risks_comparisons_male.out
nohup Rscript mp2.3_analysis_calibration_development_MSM.R > mp2.3_analysis_calibration_development_MSM.
nohup Rscript mp2.3_analysis_calibration_development_nonMSM.R > mp2.3_analysis_calibration_development_nonMSM.out
nohup Rscript mp2.3_analysis_calibration_validation_MSM.R > mp2.3_analysis_calibration_validation_MSM.
nohup Rscript mp2.3_analysis_calibration_validation_nonMSM.R > mp2.3_analysis_calibration_validation_nonMSM.out