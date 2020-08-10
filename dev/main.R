source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 32
n_h <- 10
n_s <- 10
T_thresh <- Tlength <- 115

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, check_invertability_arma)

# sample_data <- create_samples(
#   n_h = n_h, n_s = n_s, p = p, Tlength = Tlength, dim_alpha = 1,
#   percent_alpha = 0.3, range_alpha = c(0.9, 0.9),
#   ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
#   ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
#   ncores = ncores
#   )

sample_data <- prepare_corrmat_data(
  subset = 1:p,
  healthy_index_name = 'CONTROLS',

  link = "main_work/Data/Amnesia_all_AAL.mat",
  corr_matrix_name = 'corrmats',
  sick_index_name = 'TGA'

#   link = "main_work/Data/ADNI_data_AD_CN.mat",
#   corr_matrix_name = 'all.corrmats',
#   sick_index_name = 'AD'

#   link = "main_work/Data/NMDA_all_data_AAL90.mat",
#   corr_matrix_name = 'group.all',
#   sick_index_name = 'NMDA'
)

test_corr_mat(sample_data)

#iid_results <- with(sample_data$samples, estimate_loop(
#  healthy_dt = convert_corr_array_to_data_matrix_test(healthy),
#  sick_dt = convert_corr_array_to_data_matrix_test(sick)
#  ))

cov_results <- with(sample_data$samples, estimate_loop_2(
  healthy_dt = convert_corr_array_to_data_matrix_test(healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(sick)
))

plot(results$steps[[1]]$alpha, results$alpha)
abline(h = 1, v = 1)
abline(0, 1)
