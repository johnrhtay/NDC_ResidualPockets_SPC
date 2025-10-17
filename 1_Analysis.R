## ------------------------------------------------------------------------- ## 
##
## Script name: Analysis.R
## Purpose of script: Reproduce tables and figures in the manuscript
##
##    Tay JRH et al. Residual pockets as an endpoint on periodontal disease 
##    progression in supportive periodontal care. (submitted)
##
## Author: Stella Jinran Zhan
##
## ------------------------------------------------------------------------- ## 

setwd("~/NDC/John TAY")

# [0] Initialisation ---------------------------------------------------------- 
## |- Load R packages, datafile and functions ----
library(rio)
library(tidyverse)
library(logistf)
library(pROC)
library(lme4)
library(table1)
library(gtsummary)
library(rstatix)
library(predtools)

load("1_Data/Analysis datafiles/datafiles_2025-04-04.RData")
source("0_Functions.R")

## |- Parameters ----
# Outcome definitions
# progression1_2mm = Definition 1: ≥2mm PPD increase at any site (Siow et al., 2023)
# progression1a_2mm = Definition 4: At least two sites, each on a different tooth, with ≥2mm PPD increase
# progression3_3mm_2teeth = Definition 2: A ≥3mm proximal CAL loss at ≥2 teeth (Tonetti ‘standard’ definition) 
# progression4 = Definition 3: ≥1 tooth loss due to periodontal reasons
# progression5a_3mm = Definition 5: At least two sites, each on a different tooth, with ≥3mm PPD increase
confounders_vec = c("T1_Age", "Sex", "Compliance", "Staging", "Grading", "T1_%BOP", "T1_No. of teeth")
exposure_vec = "T1_No. of sites with PPD ≥5mm"
outcome_vec = c("progression1_T2_T1_2mm", "progression1a_T2_T1_2mm", "progression3_T2_T1_3mm_2teeth", "progression4_T2_T1", "progression5a_T2_T1_3mm",
                "progression1_T3_T1_2mm", "progression1a_T3_T1_2mm", "progression3_T3_T1_3mm_2teeth", "progression4_T3_T1", "progression5a_T3_T1_3mm")         
ds_analysis_red = ds_analysis %>% select(`Patient ID`, all_of(c(exposure_vec, confounders_vec, outcome_vec)), `Baseline diagnosis`, T1_n_site, contains("teeth lost")) %>%
                                   mutate(across(all_of(outcome_vec), ~ factor(., levels = c(0,1), labels = c("No", "Yes"))))              


# [1] IPCW ----
# |- Prep data for IPCW recalculate n. teeth lost based on site-level data ----
ds1_st_long = ds0_st_long %>% mutate(ID_tooth = paste0(`Patient ID`, "-", Tooth))
ID_teeth_T1 = ds1_st_long %>% filter(Time == "T1" & !is.na(PPD))

ds0_st_long_T1T2 = ds1_st_long %>% filter(Time %in% c("T1", "T2")) %>%
  filter(ID_tooth %in% ID_teeth_T1$ID_tooth) %>%
  group_by(`Patient ID`, Tooth) %>%
  reframe(n_site_PPDobserved.pertooth_T1T2 = sum(!is.na(PPD), na.rm = T)) %>% 
  mutate(observed_T1T2 = ifelse(n_site_PPDobserved.pertooth_T1T2 == 12, 1, 0))

ds0_st_long_T1T3 = ds1_st_long %>% filter(Time %in% c("T1", "T3")) %>%
  filter(ID_tooth %in% ID_teeth_T1$ID_tooth) %>%
  group_by(`Patient ID`, Tooth) %>%
  reframe(n_site_PPDobserved.pertooth_T1T3 = sum(!is.na(PPD), na.rm = T)) %>% 
  mutate(observed_T1T3 = ifelse(n_site_PPDobserved.pertooth_T1T3 == 12, 1, 0))

ds0_st_long_T1T2T3 = merge(ds0_st_long_T1T2, ds0_st_long_T1T3, by = c("Patient ID", "Tooth"))

ds_teeth_lost = ds0_st_long_T1T2T3 %>% mutate(lost_T1T2 = ifelse(observed_T1T2 == 1, 0, 1),
                                              lost_T1T3 = ifelse(observed_T1T3 == 1, 0, 1)) %>% 
  group_by(`Patient ID`) %>%
  reframe(n_tooth_lost_T1T2 = sum(lost_T1T2),
          n_tooth_lost_T1T3 = sum(lost_T1T3))

ds_analysis_red = ds_analysis %>% select(`Patient ID`, all_of(c(exposure_vec, confounders_vec, outcome_vec)), `Baseline diagnosis`, T1_n_site, contains("teeth lost")) %>%
  mutate(across(all_of(outcome_vec), ~ factor(., levels = c(0,1), labels = c("No", "Yes")))) %>% 
  right_join(ds_teeth_lost, by = "Patient ID")              

ds_summary_tooth = ds1_st_long %>% filter(Time == "T1") %>% 
  filter(ID_tooth %in% ID_teeth_T1$ID_tooth) %>%
  mutate(BOP_num = ifelse(BOP == 0, 0, 1)) %>%
  group_by(`Patient ID`, Tooth) %>%
  summarise(n_BOP = sum(BOP_num),
            n_site = sum(notNA_tooth_site_time, na.rm = T),
            `No. of sites with PPD ≥5mm` = sum(PPD >= 5, na.rm = T),
            `No. of sites with CAL ≥5mm` = sum(CAL >= 5, na.rm = T)) %>%
  mutate(`%BOP` = ifelse(is.na(n_BOP), NA, n_BOP / n_site),
         type_tooth = factor(ifelse(Tooth %in% c(17, 16, 26, 27, 37, 36, 46, 47), 1, 0),
                             levels = c(0,1), labels = c("non-molar", "molar"))) 

ds1_st_long_T1T2T3 = merge(merge(ds0_st_long_T1T2T3, ds_summary_tooth, by = c("Patient ID", "Tooth")),
                           ds_analysis_red %>% select(`Patient ID`, all_of(c(exposure_vec, confounders_vec, outcome_vec)), contains("lost")), by = "Patient ID")


# |- IPCW model ----
ipcw_data = ds1_st_long_T1T2T3 
IPCW_var = c("No. of sites with PPD ≥5mm", "type_tooth", "T1_Age", "Staging", "Grading", "Compliance", "Sex",
             "No. of sites with CAL ≥5mm", "%BOP", "T1_No. of teeth")

table1(~ factor(observed_T1T2) + factor(observed_T1T3), data = ipcw_data)

## |- T1 to T2 ----
# Covariates for the IPCW model (T1 to T2): 
# - number of sites with PPDs ≥5mm
# - tooth type (molar versus non-molar)
# - age, sex
# - staging, grading, compliance, 
# - number of sites with CAL ≥5 mm, 
# - total number of teeth at T1 
out1_T1T2 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                     outcome_var_IPCW = "observed_T1T2", teeth_lost_period = "teeth lost (T1-T2)", IPCW_var_m = IPCW_var[-9], outcome_var = outcome_vec[1])
out1a_T1T2 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                      outcome_var_IPCW = "observed_T1T2", teeth_lost_period = "teeth lost (T1-T2)",IPCW_var_m = IPCW_var[-9], outcome_var = outcome_vec[2])
out3_T1T2 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                     outcome_var_IPCW = "observed_T1T2", teeth_lost_period = "teeth lost (T1-T2)", IPCW_var_m = IPCW_var[-9], outcome_var = outcome_vec[3])
out4_T1T2 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                     outcome_var_IPCW = "observed_T1T2", teeth_lost_period = "teeth lost (T1-T2)", IPCW_var_m = IPCW_var[-9], outcome_var = outcome_vec[4])
out5a_T1T2 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                      outcome_var_IPCW = "observed_T1T2", teeth_lost_period = "teeth lost (T1-T2)", IPCW_var_m = IPCW_var[-9], outcome_var = outcome_vec[5])

## |- T1 to T3 ----
# Covariates for the IPCW model (T1 to T3), reduced model due to no convergence: 
# - number of sites with PPDs ≥5mm
# - tooth type (molar versus non-molar)
# - age, sex
# - staging, compliance
out1_T1T3 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                     outcome_var_IPCW = "observed_T1T3", teeth_lost_period = "teeth lost (T1-T3)", IPCW_var_m = IPCW_var[c(1:4, 6:7)], outcome_var = outcome_vec[6])
out1a_T1T3 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                      outcome_var_IPCW = "observed_T1T3", teeth_lost_period = "teeth lost (T1-T3)", IPCW_var_m = IPCW_var[c(1:4, 6:7)], outcome_var = outcome_vec[7])
out3_T1T3 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                     outcome_var_IPCW = "observed_T1T3", teeth_lost_period = "teeth lost (T1-T3)", IPCW_var_m = IPCW_var[c(1:4, 6:7)], outcome_var = outcome_vec[8])
out4_T1T3 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                     outcome_var_IPCW = "observed_T1T3", teeth_lost_period = "teeth lost (T1-T3)", IPCW_var_m = IPCW_var[c(1:4, 6:7)], outcome_var = outcome_vec[9])
out5a_T1T3 = fct_IPCW(ds_analysis_pt_level=ds_analysis_red, ipcw_data=ipcw_data, exposure_vec=exposure_vec, confounders_vec=confounders_vec,
                      outcome_var_IPCW = "observed_T1T3", teeth_lost_period = "teeth lost (T1-T3)", IPCW_var_m = IPCW_var[c(1:4, 6:7)], outcome_var = outcome_vec[10])

ds_MV_list = list(out1_T1T2, out1a_T1T2, out3_T1T2, out4_T1T2, out5a_T1T2,
                  out1_T1T3, out1a_T1T3, out3_T1T3, out4_T1T3, out5a_T1T3)


# [2] ROC analysis ----
# Models:
# - main-effect models for 3, 4, and 5a
# - with and without interactions for 1 and 1a

## |- ROC curves for main-effect models ----
roc_progression_model_list = tb_roc_list = list()
for (i in 1:10){
  roc_progression_model_list[[i]] = fct_ROC(outcome_var = outcome_vec[i], ds_ROC = ds_MV_list[[i]]$ds_final_out)
  tb_roc_list[[i]] = roc_progression_model_list[[i]]$tb_roc
  tiff(paste0("Plots/ROC_", outcome_vec[i],"_", today(), ".tiff"), width = 8, height = 7, units="in", res=300)
  plot(roc_progression_model_list[[i]]$roc_progression_MV, grid=T, print.thres = T, print.auc = T, col="blue", main = outcome_vec[i])
  dev.off()
}

# [3] Optimal cut-off for exposure ----
roc_cutoff_list = list()
tb_cutoff_list = c()
for (i in 1:10){
  roc_cutoff_list[[i]] = fct_ROC_cutoff(outcome_var = outcome_vec[i])
  tb_cutoff = roc_cutoff_list[[i]]$tb_cutoff
  tb_cutoff_list = rbind(tb_cutoff_list, tb_cutoff)
}


# [4] Calibration plots ----
for (i in 1:10){
  ds_cal = ds_MV_list[[i]]$ds_final_out %>% mutate(outcome = ifelse(get(outcome_vec[i]) == "Yes", 1, 0))

  tiff(paste0("Plots/CalibrationPlot_IPCW_", outcome_vec[i],"_", today(), ".tiff"), width = 8, height = 7, units="in", res=300)
  fct_calibration(model_name = outcome_vec[i], pred_var = "pred_MV", ds_cal = ds_cal)
  dev.off()
}
