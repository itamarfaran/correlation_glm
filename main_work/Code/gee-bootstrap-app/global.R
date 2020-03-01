sidebar_options <- list(
  'None' = 'none',
  '1' = '1',
  'P' = 'p',
  'N' = 'n',
  'NH' = 'nh',
  'ND' = 'nd',
  'P/N' = 'p_n_ratio',
  'Length of Times Serie' = 'Tlength',
  'Percentage' = 'sick_obs_percentage',
  'Actual SD' = 'actual_sd',
  'GEE Estimated SD' = 'gee_sd',
  'GEE(old) Estimated SD' = 'gee_old_sd',
  'MLE Estimated SD' = 'mle_sd',
  'Ratio' = 'actual_est_ratio',
  'Type I Error' = 'type1error'
)

sidebar_options_reverse <- list()
for(i in 1:length(sidebar_options)){
  sidebar_options_reverse[[i]] <- names(sidebar_options)[i]
  names(sidebar_options_reverse)[i] <- sidebar_options[[i]]
}
