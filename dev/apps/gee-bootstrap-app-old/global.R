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
  'Actual SD' = 'actual_sd',
  'GEE Estimated SD' = 'gee_sd',
  'GEE(old) Estimated SD' = 'gee_old_sd',
  'MLE Estimated SD' = 'mle_sd',
  'Ratio' = 'actual_est_ratio',
  'Type I Error' = 'type1error'
)

sidebar_options_reverse <- list()
for(i in seq_along(sidebar_options)){
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
