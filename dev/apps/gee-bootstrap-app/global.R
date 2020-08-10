library(ggplot2)

sidebar_options <- list(
  'None' = 'none',
  '1' = '1',
  'P' = 'p',
  'M' = 'm',
  'N' = 'n',
  'NH' = 'nh',
  'ND' = 'nd',
  'AR Coefficient' = 'ar',
  'P/N' = 'p_n_ratio',
  'M/N' = 'm_n_ratio',
  'Length of Times Serie' = 'Tlength',
  'Percentage of Sick Samples' = 'sick_obs_percentage',
  'None Null Alphas' = 'percent_alpha',
  'Minimum Range of Alpha' = 'min_alpha',
  'Real Alpha' = 'alpha',
  'Alpha Estimate' = 'alpha_est_mean',
  'Alpha Estimate SD' = 'alpha_est_emp_sd',
  'Alpha Estimate Lower' = 'alpha_est_lower',
  'Alpha Estimate Upper' = 'alpha_est_upper',
  'GEE Estimated SD' = 'gee_sd_mean',
  'SD of GEE Estimated SD' = 'gee_sd_emp_sd',
  'GEE Estimated SD Lower' = 'gee_sd_lower',
  'GEE Estimated SD Upper' = 'gee_sd_upper')

sidebar_options_reverse <- list()
for(i in 1:length(sidebar_options)){
  sidebar_options_reverse[[i]] <- names(sidebar_options)[i]
  names(sidebar_options_reverse)[i] <- sidebar_options[[i]]
}

scales_x <- list(
  'Continuous' = scale_x_continuous,
  'Log10' = scale_x_log10,
  'Reverse' = scale_x_reverse,
  'Sqrt' = scale_x_sqrt
)

scales_y <- list(
  'Continuous' = scale_y_continuous,
  'Log10' = scale_y_log10,
  'Reverse' = scale_y_reverse,
  'Sqrt' = scale_y_sqrt
)
