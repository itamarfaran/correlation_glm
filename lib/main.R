source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 22
n_h <- 10
n_s <- 10
T_thresh <- Tlength <- 115

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, check_invertability_arma)

sample_data <- create_samples(
  n_h = n_h, n_s = n_s, p = p, Tlength = Tlength, dim_alpha = 1,
  percent_alpha = 0.3, range_alpha = c(0.9, 0.9),
  ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
  ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
  ncores = ncores
  )

# sample_data <- prepare_corrmat_data(
#   subset = 1:p,
#   healthy_index_name = 'CONTROLS',

#   link = "lib/data/Amnesia_all_AAL.mat",
#   corr_matrix_name = 'corrmats',
#   sick_index_name = 'TGA'

#   link = "lib/data/ADNI_data_AD_CN.mat",
#   corr_matrix_name = 'all.corrmats',
#   sick_index_name = 'AD'

#   link = "lib/data/NMDA_all_data_AAL90.mat",
#   corr_matrix_name = 'group.all',
#   sick_index_name = 'NMDA'
# )

test_corr_mat(sample_data)

results <- with(sample_data$samples, estimate_alpha(healthy_dt = healthy, sick_dt = sick))

gee_var <- with(sample_data$samples, compute_gee_variance(
  healthy_dt = healthy, sick_dt = sick, cov_obj = results
  ))

results_jacknife <- estimate_alpha_jacknife(
  healthy_dt = sample_data$samples$healthy, sick_dt = sample_data$samples$sick,
  linkFun = linkFun, jack_healthy = TRUE, return_gee = TRUE, ncores = ncores)

alpha_jk_estimate <- with(results_jacknife, colMeans(alpha))
alpha_jk_variance <- with(results_jacknife, var(alpha)*(nrow(alpha) - 1)^2/nrow(alpha)) 
alpha_jk_sd_estimate <- t(apply(
  X = results_jacknife$gee_var, MARGIN = 1,
  FUN = function(x) sqrt(diag(vector2triangle(x, diag = TRUE)))
  ))

# (n_h + n_s - 1)*sqrt_diag(gee_var) - colMeans(alpha_jk_sd_estimate)


out <- data.table(
  alpha = as.vector(sample_data$alpha),
  null = as.vector(sample_data$alpha == linkFun$NULL_VAL),
  gee_estimate = as.vector(results$alpha),
  gee_sd = sqrt_diag(gee_var),
  jk_estimate = alpha_jk_estimate,
  jk_sd = sqrt_diag(alpha_jk_variance)
)

out[,.(
  gee_estimate_sd = sd(gee_estimate),
  gee_sd_estimate = mean(gee_sd),
  jk_estimate_sd = sd(jk_estimate),
  jk_sd_estimate = sd(jk_estimate)
  ),
  by = null]

out[,.(
  alpha, null,
  z_gee = (gee_estimate - linkFun$NULL_VAL)/gee_sd,
  z_jk = (jk_estimate - linkFun$NULL_VAL)/jk_sd
)][,.(
  alpha, null,
  p_gee = 2*pnorm(abs(z_gee), lower.tail = F),
  p_jk = 2*pnorm(abs(z_jk), lower.tail = F)
)][,.(
  gee_sig_01 = mean(p_gee < 0.1),
  gee_sig_005 = mean(p_gee < 0.05),
  gee_sig_001 = mean(p_jk < 0.01),
  jk_sig_01 = mean(p_jk < 0.1),
  jk_sig_005 = mean(p_jk < 0.05),
  gee_sig_001 = mean(p_gee < 0.01)
), keyby = null]



mean_sqrt_diag(alpha_jk_variance)
sd(results$alpha)
mean_sqrt_diag(gee_var)

# save.image('main_work/data/enviroments/full_run_jacknife_simulated_data_multiplicative_link.RData')
