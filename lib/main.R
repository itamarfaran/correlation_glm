source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 200
n_h <- 10
n_s <- 10
T_thresh <- Tlength <- 115

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, check_invertability_arma)

#sample_data <- create_samples(
#  n_h = n_h, n_s = n_s, p = p, Tlength = Tlength, dim_alpha = 1,
#  percent_alpha = 0.3, range_alpha = c(0.9, 0.9),
#  ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
#  ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
#  ncores = ncores
#  )

 sample_data <- prepare_corrmat_data(
   subset = 1:p,
   healthy_index_name = 'CONTROLS',

#  link = "lib/data/Amnesia_all_AAL.mat",
   link = "lib/data/Amnesia_all_AICHA.mat",
   corr_matrix_name = 'corrmats',
   sick_index_name = 'TGA'

#   link = "lib/data/ADNI_data_AD_CN.mat",
#   corr_matrix_name = 'all.corrmats',
#   sick_index_name = 'AD'

#   link = "lib/data/NMDA_all_data_AAL90.mat",
#   corr_matrix_name = 'group.all',
#   sick_index_name = 'NMDA'
 )

test_corr_mat(sample_data)

results <- with(sample_data$samples, estimate_alpha(healthy_dt = healthy, sick_dt = sick, early_stop = TRUE))

gee_var <- with(sample_data$samples, compute_gee_variance(
  healthy_dt = healthy, sick_dt = sick, cov_obj = results
  ))

# results_jacknife <- estimate_alpha_jacknife(
#   healthy_dt = sample_data$samples$healthy, sick_dt = sample_data$samples$sick,
#   linkFun = linkFun, jack_healthy = TRUE, return_gee = F, ncores = ncores)
# jacknife_inference <- infer_jacknife(results_jacknife)

save.image('main_work/data/enviroments/full_run_jacknife_simulated_data_multiplicative_link.RData')

