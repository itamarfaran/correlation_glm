sidebar_options <- list(
  'None' = 'none',
  'P' = 'p',
  'N' = 'n',
  'Percentage' = 'sick_obs_percentage',
  'Actual SD' = 'actual_sd',
  'GEE Estimated SD' = 'mean_est_sd',
  'Ratio' = 'actual_est_ratio'
)

sidebar_options_reverse <- list()
for(i in 1:length(sidebar_options)){
  sidebar_options_reverse[[i]] <- names(sidebar_options)[i]
  names(sidebar_options_reverse)[i] <- sidebar_options[[i]]
}
 