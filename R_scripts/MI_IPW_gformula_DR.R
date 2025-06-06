if (!require('tidyverse')) install.packages('tidyverse')
library(tidyverse)
if (!require('mice')) install.packages('mice')
library(mice)
if (!require('modelr')) install.packages('modelr')
library(modelr)

data_dir <- './data/'
output_dir <- './output/'

version <- '250430'

df_crude <- read.csv(paste0(data_dir, version, '_covariates_joined.csv'))

df <- df_crude %>%
  mutate(vasopressin = ifelse(first_use_drug == 'vasopressin', 1, 0)) %>% 
  select(-c(first_use_drug, time_zero, primary_diagnosis, icu_stay_id,
            past_rrt, past_vent,
            period_intime_to_timezero, aki_stage, hospital_admission_id,
            aki_onset, aki_worse, egfr_worse, make_outcome))

imp.meth <- map_chr(df, ~ ifelse(any(is.na(.x)), "norm", ""))
print(imp.meth)

n_imputations <- 100
list_df_imp <- mice(df, m=n_imputations, maxit=50, meth=imp.meth, seed=813)

set.seed(813)

for (i in 1:n_imputations) {
  df_imp <- mice::complete(list_df_imp, i)
  assign(paste0("df_", i), df_imp)
  rm(df_imp)
  save(list = paste0("df_", i), file = paste0(data_dir, "mi_data/imp", i, ".RData"))
}

ipw_dr_analysis <- function(df) {

  # Estimation of propensity score
  treat_model <- glm(vasopressin ~ female + age + I(age^2) + cci + period_from_first_norad +
                       apache2 + I(apache2^2) + sofa + lactate + creatinine +
                       on_vent  + on_rrt  + standardized_norad_rate + I(standardized_norad_rate^2) +
                       pf_ratio + infected_nervous_system + infected_cardiovascular + infected_respiratory +
                       infected_abdomen + infected_urinary_tract + infected_soft_tissue + infected_other +
                       hr_pre + mbp_pre + I(mbp_pre^2) + bt_pre,
                     data = df,
                     family = binomial())
  
  df <- df %>% mutate(
    prob_exposure = predict(treat_model, df, type = 'response'),
    crude_w = ifelse(vasopressin == 1, 1/prob_exposure, 1/(1-prob_exposure))
  )
  
  w_percentile_99.5 <- df %>%
    pull(crude_w) %>%
    quantile(., 0.995)
  
  df <- df %>% mutate(
    w = pmin(crude_w, w_percentile_99.5)
  )
  
  # Marginal Structural Model
  msm_in_hospital_death <- lm(in_hospital_death ~ vasopressin,
                              data = df,
                              weights = w)
  
  rd_in_hospital_death_ipw <- coef(msm_in_hospital_death)[2]
  
  # Doubly Robust Estimator
  df <- df %>% 
    mutate(r = ifelse(vasopressin == 1, w, -w))
  
  df$interv <- -1
  
  df_0 <- df %>% 
    mutate(
      interv = 0,
      vasopressin = 0,
      in_hospital_death = NA,
      crude_w = 1/(1-prob_exposure)
    )
  
  df_0_w_percentile_99.5 <- df_0 %>%
    pull(crude_w) %>%
    quantile(., 0.995)
  
  df_0 <- df_0 %>% 
    mutate(
      w = pmin(crude_w, df_0_w_percentile_99.5),
      r = -w
    )
  
  df_1 <- df %>% 
    mutate(
      interv = 1,
      vasopressin = 1,
      in_hospital_death = NA,
      crude_w = 1/prob_exposure
    )
  
  df_1_w_percentile_99.5 <- df_1 %>%
    pull(crude_w) %>%
    quantile(., 0.995)
  
  df_1 <- df_1 %>% 
    mutate(
      w = pmin(crude_w, df_1_w_percentile_99.5),
      r = w
    )
  
  df_onesample <- rbind(df, df_0, df_1)
  
  dr_model_in_hospital_death <- glm(in_hospital_death ~ vasopressin*standardized_norad_rate + vasopressin*I(standardized_norad_rate^2) +
                                      vasopressin*mbp_pre + vasopressin*I(mbp_pre^2) +
                                      female + age + I(age^2) + cci + period_from_first_norad +
                                      apache2 + I(apache2^2) + sofa + lactate + creatinine +
                                      on_vent  + on_rrt +
                                      pf_ratio + infected_nervous_system + infected_cardiovascular + infected_respiratory +
                                      infected_abdomen + infected_urinary_tract + infected_soft_tissue + infected_other +
                                      hr_pre + bt_pre +
                                      r,
                                    data=df_onesample,
                                    family=binomial())
  
  df_onesample <- df_onesample %>% 
    mutate(prob = predict(dr_model_in_hospital_death, df_onesample, type = 'response'))
  
  untreated_countfac <- df_onesample %>% filter(interv==0) %>% pull(prob) %>% mean
  treated_countfac <- df_onesample %>% filter(interv==1) %>% pull(prob) %>% mean
  
  rd_in_hospital_death_dr<- treated_countfac - untreated_countfac
  
  return(c(rd_in_hospital_death_ipw, rd_in_hospital_death_dr))
}

g_formula_analysis <- function(df) {
  df$interv <- -1
  
  df_0 <- df %>% 
    mutate(
      interv = 0,
      vasopressin = 0,
      in_hospital_death = NA
    )
  
  df_1 <- df %>% 
    mutate(
      interv = 1,
      vasopressin = 1,
      in_hospital_death = NA
    )

  df_onesample <- rbind(df, df_0, df_1)
  
  g_outcome_model_in_hospital_death <- glm(in_hospital_death ~ vasopressin*standardized_norad_rate + vasopressin*I(standardized_norad_rate^2) +
                                             vasopressin*mbp_pre + vasopressin*I(mbp_pre^2) +
                                             female + age + I(age^2) + cci + period_from_first_norad +
                                             apache2 + I(apache2^2) + sofa + lactate + creatinine +
                                             on_vent  + on_rrt +
                                             pf_ratio + infected_nervous_system + infected_cardiovascular + infected_respiratory +
                                             infected_abdomen + infected_urinary_tract + infected_soft_tissue + infected_other +
                                             hr_pre + bt_pre,
                         data=df_onesample,
                         family=binomial())
  
  df_onesample <- df_onesample %>% 
    mutate(prob = predict(g_outcome_model_in_hospital_death, df_onesample, type = 'response'))
  
  untreated_countfac <- df_onesample %>% filter(interv==0) %>% pull(prob) %>% mean
  treated_countfac <- df_onesample %>% filter(interv==1) %>% pull(prob) %>% mean
  
  rd_in_hospital_death_gformula <- treated_countfac - untreated_countfac
  return(rd_in_hospital_death_gformula)
}

bootstrap_in_hospital_death <- function(df, n_iter){
  rd_in_hospital_ipw <- as.vector(NULL)
  rd_in_hospital_gformula <- as.vector(NULL)
  rd_in_hospital_dr <- as.vector(NULL)

  for (i in 1:n_iter){
    resample_data <- resample_bootstrap(df)
    resample_df <- as.data.frame(resample_data)
    
    rd_ipw_dr <- ipw_dr_analysis(resample_df)
    rd_in_hospital_ipw[i] <- rd_ipw_dr[1]
    rd_in_hospital_gformula[i] <- g_formula_analysis(resample_df)
    rd_in_hospital_dr[i] <- rd_ipw_dr[2]
  }
  
  mean_ipw <- rd_in_hospital_ipw %>% mean
  var_ipw <- rd_in_hospital_ipw %>% var
  mean_gformula <- rd_in_hospital_gformula %>% mean
  var_gformula <- rd_in_hospital_gformula %>% var
  mean_dr <- rd_in_hospital_dr %>% mean
  var_dr <- rd_in_hospital_dr %>% var
  return(c(mean_ipw, var_ipw, mean_gformula, var_gformula, mean_dr, var_dr))
}

mean_imp_ipw <- as.vector(NULL)
var_imp_ipw <- as.vector(NULL)
mean_imp_gformula <- as.vector(NULL)
var_imp_gformula <- as.vector(NULL)
mean_imp_dr <- as.vector(NULL)
var_imp_dr <- as.vector(NULL)

n_iter <- 1000

# Starting time
time1 <- Sys.time()

for (i in 1:n_imputations) {
  load_path <- paste0(data_dir, "mi_data/imp", i, ".RData")
  load(load_path)
  
  df_name <- paste0("df_", i)
  df_imp <- get(df_name)
  
  df_imp <- df_imp %>% 
    mutate(standardized_norad_rate = norad_rate_pre / body_weight_pre) %>% 
    select(-c(norad_rate_pre, body_weight_pre))
  
  results <- bootstrap_in_hospital_death(df_imp, n_iter)
  
  mean_imp_ipw[i] <- results[1]
  var_imp_ipw[i] <- results[2]
  mean_imp_gformula[i] <- results[3]
  var_imp_gformula[i] <- results[4]
  mean_imp_dr[i] <- results[5]
  var_imp_dr[i] <- results[6]
  
  rm(df_imp)
}

time2 <- Sys.time()
# Measure time to complete analysis
difftime(time2, time1, units='mins')

rubins_rule <- function(mean_imp, var_imp, n_imputations){
  ATE <- mean(mean_imp)
  WM <- mean(var_imp)
  sigma <- as.vector(NULL)
  for (i in 1:n_imputations){
    sigma[i] <- (mean_imp[[i]] - ATE)^2
  }
  BM <- sum(sigma)/(n_imputations - 1)
  TM <- WM + (1+(1/n_imputations))*BM
  se <- sqrt(TM)
  ll <- ATE - qnorm(0.975)*se
  ul <- ATE + qnorm(0.975)*se
  results <- data.frame(ATE, ll, ul, se)
  return(results)
}

ipw_results <- rubins_rule(mean_imp_ipw, var_imp_ipw, n_imputations)
gformula_results <- rubins_rule(mean_imp_gformula, var_imp_gformula, n_imputations)
dr_results <- rubins_rule(mean_imp_dr, var_imp_dr, n_imputations)

if (!require('forestplot')) install.packages('forestplot')
library(forestplot)

results_fp <- data.frame(
  methods = c('IP weighting', 'g-formula', 'Doubly Robust Estimator'),
  sample_size = c(nrow(df), nrow(df), nrow(df)),
  estimates = c(ipw_results$ATE, gformula_results$ATE, dr_results$ATE),
  lwr = c(ipw_results$ll, gformula_results$ll, dr_results$ll),
  upr = c(ipw_results$ul, gformula_results$ul, dr_results$ul),
  rd = c(
    paste0(format(100 * round(ipw_results$ATE, 3), nsmall = 1), "% (",
           format(100 * round(ipw_results$ll, 3), nsmall = 1), "% to ",
           format(100 * round(ipw_results$ul, 3), nsmall = 1), "%)"),
    paste0(format(100 * round(gformula_results$ATE, 3), nsmall = 1), "% (",
           format(100 * round(gformula_results$ll, 3), nsmall = 1), "% to ",
           format(100 * round(gformula_results$ul, 3), nsmall = 1), "%)"),
    paste0(format(100 * round(dr_results$ATE, 3), nsmall = 1), "% (",
           format(100 * round(dr_results$ll, 3), nsmall = 1), "% to ",
           format(100 * round(dr_results$ul, 3), nsmall = 1), "%)")
  )
)

forest_plot <- results_fp %>% 
  forestplot(
    mean = estimates,
    lower = lwr,
    upper = upr,
    labeltext = c(methods, sample_size, rd),
    xlog = FALSE,
    sizes = 0.75,
    zero = 0,
    title = 'Effect of vasopressin as a second-line agents on in-hospital mortality',
    xticks = seq(-0.1, 0.1, by = 0.05) 
  ) %>% 
  fp_set_style(
    box = 'black',
    line = 'black',
    align = 'lcc',
    txt_gp = fpTxtGp(ticks = gpar(cex = 1))
  ) %>% 
  fp_add_header(
    methods = 'Adjustment methods',
    sample_size = 'N',
    rd = 'Risk difference (95% CI)'
  ) %>% 
  fp_decorate_graph(graph.pos = 3) %>% 
  fp_set_zebra_style('#EFEFEF')
forest_plot

output_file <- paste0(output_dir, 'forest_plot.png')
png(output_file, width = 10, height = 6, units = 'in', res = 300)
print(forest_plot)
dev.off()
