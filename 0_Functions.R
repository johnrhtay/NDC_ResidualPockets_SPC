## ------------------------------------------------------------------------- ## 
##
## Script name: Functions.R
## Purpose of script: Functions to perform the analyses in the manuscript
##
##    Tay JRH et al. Residual pockets as an endpoint on periodontal disease 
##    progression in supportive periodontal care. (submitted)
##
## Author: Stella Jinran Zhan
##
## ------------------------------------------------------------------------- ## 

# fct_IPCW: to obtain IPCW-weighted logistic regression results
# Parameters: 
#   @ds_analysis_pt_level = analysis dataset (patient-level)
#   @ipcw_data = dataset with info on observed tooth (tooth-level)
#   @outcome_var_IPCW = outcome variable for IPCW (i.e observed_T1T2 or observed_T1T3)
#   @IPCW_var_m = covariates to be included in the IPCW model for each outcome
#   @outcome_var = outcome variable for the final regression analysis on progression
#   @exposure_vec = T1_No. of sites with PPD ≥5mm
#   @confounders_vec = covariates for the final regression analysis on progression
#   @teeth_lost_period = "teeth lost (T1-T2)" or "teeth lost (T1-T3)"
#   @trim = indicator variable for trimming extreme values in IPCW
fct_IPCW = function(ds_analysis_pt_level, ipcw_data, outcome_var_IPCW, IPCW_var_m, 
                    outcome_var, exposure_vec, confounders_vec, teeth_lost_period = "teeth lost (T1-T2)", trim = F){
  
  # [1] Set up IPCW ----
  ## Calculate tooth-level weights from IPCW
  m_formula = paste0(outcome_var_IPCW, " ~ `", paste0(IPCW_var_m, collapse = "` + `"), "` + (1 | `Patient ID`)")
  IPCW_model = glmer(m_formula, data = ipcw_data, family = binomial)
  
  ## ipcw_data = ipcw_data %>% filter(!is.na(`%BOP`))
  ipcw_data$IPCW = 1/predict(IPCW_model, newdata = ipcw_data, type = "response")
  ## Inspect weights
  print(summary(ipcw_data$IPCW))
  hist(ipcw_data$IPCW)
  # boxplot(ipcw_data$IPCW)
  
  ## Testing trimming extreme values
  if(trim == T){
    trim_thres = quantile(ipcw_data$IPCW, 0.99)
    ipcw_data$IPCW_trimmed = pmin(ipcw_data$IPCW, trim_thres)
    print(summary(ipcw_data$IPCW_trimmed))
    ipcw_data$IPCW = ipcw_data$IPCW_trimmed
  }
  
  ## Calculate patient-level weights from the tooth-level using only the observed data
  ipcw_data.cc = ipcw_data %>% filter(get(outcome_var_IPCW) == 1)
  IPCW_average_pt = ipcw_data.cc %>% group_by(`Patient ID`) %>%
    summarise(av.weight_pt = mean(IPCW, na.rm = TRUE)) %>% ungroup()
  
  # [2] Create outputs ----
  ## Descriptive
  tbl_summary = ds_analysis_red %>% select(all_of(c(exposure_vec, confounders_vec)), outcome_var, 
                                           `Baseline diagnosis`, T1_n_site, contains(teeth_lost_period)) |>
    tbl_summary(by = outcome_var,
                digits = list(all_categorical() ~ c(0, 1),
                              all_continuous() ~ 0))
  ## Univariable
  UNIreg = ds_analysis %>% select(all_of(c(exposure_vec, confounders_vec, outcome_var))) %>% 
    tbl_uvregression(method = glm, y = outcome_var, method.args = list(family = binomial),  
                     exponentiate = TRUE, hide_n = TRUE,
                     pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
    # add_nevent() %>% 
    bold_p() %>% bold_labels() %>% italicize_levels() %>% 
    modify_table_styling(column = estimate, rows = !is.na(estimate),
                         cols_merge_pattern = "{estimate} ({conf.low})") %>%
    modify_header(estimate ~ "**OR (95% CI)**") 
  
  ## Run IPCW-weighted logistic regression
  ds_final_weight = merge(ds_analysis_pt_level, IPCW_average_pt, by = "Patient ID")
  MV_formula = as.formula(c(exposure_vec, confounders_vec) %>% str_c(collapse = "` + `") %>% str_c(outcome_var, " ~ `", ., "`") )
  
  MVreg = glm(MV_formula, family = "binomial", data = ds_final_weight)
  MVreg_weighted = glm(MV_formula, family = "quasibinomial", data = ds_final_weight, weights = av.weight_pt)
  
  MVreg_tbl = tbl_regression(MVreg, exponentiate = TRUE, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
    bold_p() %>% bold_labels() %>% italicize_levels() %>% add_nevent() %>%
    modify_table_styling(column = estimate, rows = !is.na(estimate),
                         cols_merge_pattern = "{estimate} ({conf.low})") %>%
    modify_header(estimate ~ "**OR (95% CI)**")
  
  MVreg_weighted_tbl = tbl_regression(MVreg_weighted, exponentiate = TRUE, pvalue_fun = ~ style_pvalue(.x, digits = 3)) %>%
    bold_p() %>% bold_labels() %>% italicize_levels() %>% add_nevent() %>%
    modify_table_styling(column = estimate, rows = !is.na(estimate),
                         cols_merge_pattern = "{estimate} ({conf.low})") %>%
    modify_header(estimate ~ "**OR (95% CI)**")
  
  theme_gtsummary_compact()
  MV_out = tbl_merge(list(UNIreg, MVreg_tbl, MVreg_weighted_tbl),
                     tab_spanner = c(paste0("**Univariable:", outcome_var, "**"),
                                     paste0("**Unweighted:", outcome_var, "**"),
                                     paste0("**Weighted:", outcome_var, "**")))
  
  model_fit = tibble(Parameter = c(exposure_vec, confounders_vec, "Hosmer-Lemeshow Test"),
                     Value = c(car::vif(MVreg_weighted), ResourceSelection::hoslem.test(MVreg_weighted$y, fitted(MVreg_weighted))$p.value))
  
  ## Calibration plot
  library(rms)
  fit_lrm = lrm(MV_formula, x = TRUE, y = TRUE, data = ds_final_weight, weights = av.weight_pt)
  
  ## final dataset with predicted values from final weighted model
  ds_final_out = ds_final_weight %>% mutate(pred_MV = predict(MVreg_weighted, type = "response"))
  
  return(list(tbl_summary=tbl_summary, MV_out=MV_out, model_fit=model_fit, fit_lrm=fit_lrm, ds_final_out=ds_final_out))
}

# fct_ROC: ROC analysis
# Parameters: 
#   @outcome_var = outcome variable
#   @exposure_var = exposure variable
#   @ds_ROC = analysis dataset with predicted prob from the final model
fct_ROC = function(outcome_var, exposure_var = exposure_vec, ds_ROC){
  # outcome_var = outcome_vec[1]
  
  # roc(actual binary outcome vs predicted prob from the model)
  roc_progression_MV = roc(response = ds_ROC[[outcome_var]], predictor = ds_ROC[["pred_MV"]], direction = "<", ci=TRUE, auc=TRUE)

  tb_roc = tibble(outcome = outcome_var, model = "MV", 
                  auc = auc(roc_progression_MV),
                  auc_CI = paste0("(", round(ci.auc(roc_progression_MV)[1], 3), ", ", round(ci.auc(roc_progression_MV)[3], 3), ")"))
  
  return(list(tb_roc=tb_roc, roc_progression_MV=roc_progression_MV))
}

# fct_ROC_cutoff: optimal cutoff for exposure variable
# Parameters: 
#   @outcome_var = outcome variable
#   @exposure_var = exposure variable
fct_ROC_cutoff = function(outcome_var, exposure_var = exposure_vec){
  
  roc_progression = roc(ds_analysis[[outcome_var]], ds_analysis[[exposure_var]], direction = "<", ci=T, auc=T)
  optimal_cutoff_closest = coords(roc_progression, "best", ret = c("threshold", "sensitivity", "specificity", "npv", "ppv"), best.method="closest.topleft")
  optimal_cutoff_youden = coords(roc_progression, "best", ret = c("threshold", "sensitivity", "specificity", "npv", "ppv"), best.method="youden")
  tb_cutoff = c()
  for (cutoff_n in 1:10){
    tb_cutoff_n = coords(roc_progression, cutoff_n, ret = c("threshold", "sensitivity", "specificity", "npv", "ppv"))
    tb_cutoff = rbind(tb_cutoff, tb_cutoff_n)
  }
  tb_cutoff = rbind(tb_cutoff, optimal_cutoff_closest, optimal_cutoff_youden)
  tb_cutoff$type = c(rep("fixed value", 10), "closest.topleft", "youden")
  tb_cutoff$youden = tb_cutoff$sensitivity + tb_cutoff$specificity - 1
  tb_cutoff$auc = auc(roc_progression)
  tb_cutoff$auc.cil = ci.auc(roc_progression)[1]
  tb_cutoff$auc.ciu = ci.auc(roc_progression)[3]
  tb_cutoff = tb_cutoff %>% mutate(outcome = outcome_var, 
                                   auc_CI = paste0("(", round(auc.cil, 3), ", ", round(auc.ciu, 3), ")")) %>%
    select(outcome, type, threshold:ppv, youden, auc, auc_CI)
  return(list(tb_cutoff=tb_cutoff, roc_progression=roc_progression))
}

# fct_calibration: calibration plot
# Reference: https://darrendahly.github.io/post/homr/
# Parameters: 
#   @model_name = plot title (e.g. outcome variable )
#   @pred_var = column containing the predicted prob from the final model
#   @ds_cal = analysis dataset with predicted prob from the final model
fct_calibration <- function(model_name, pred_var, ds_cal){
  
  require(tidyverse)
  require(viridis)
  require(gridExtra)
  
  # The calibration plot        
  g1 <- mutate(ds_cal, bin = ntile(get(pred_var), 10)) %>% 
    # Bin prediction into 10ths
    group_by(bin) %>%
    mutate(n = n(), # Get ests and CIs
           bin_pred = mean(get(pred_var)), 
           bin_prob = mean(as.numeric(outcome)), 
           se = sqrt((bin_prob * (1 - bin_prob)) / n), 
           ul = bin_prob + 1.96 * se, 
           ll = bin_prob - 1.96 * se) %>%
    ungroup() %>%
    ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
    geom_pointrange(size = 0.5, color = "black") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    geom_abline() + # 45 degree line indicating perfect calibration
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", 
                color = "black") + 
    # straight line fit through estimates
    geom_smooth(aes(x = get(pred_var), y = as.numeric(outcome)), 
                color = "red", se = FALSE, method = "loess") + 
    # loess fit through estimates
    xlab("") +
    ylab("Observed Probability") +
    theme_minimal() +
    ggtitle(model_name)
  
  # The distribution plot        
  g2 <- ggplot(ds_cal, aes(x = get(pred_var))) +
    geom_histogram(fill = "black", bins = 200) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    xlab("Predicted Probability") +
    ylab("") +
    theme_minimal() +
    scale_y_continuous(breaks = c(0, 40)) +
    theme(panel.grid.minor = element_blank())
  
  # Combine them    
  g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
  grid::grid.newpage()
  grid::grid.draw(g)
  return(g[[3]])
  
}
