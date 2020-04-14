source("main_work/code/01_general_functions.R")
source("main_work/code/02_simulation_functions.R")
source("main_work/code/03_estimation_functions.R")
source("main_work/code/04_inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 20
nH <- 30
nS <- 30
T_thresh <- Tlength <- 115

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, check_invertability_arma)

sampleData <- create_samples(
  nH = nH, nS = nS, p = p, Tlength = Tlength, dim_alpha = 1,
  percent_alpha = 0.3, range_alpha = c(1, 1),
  ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
  ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
  ncores = ncores
  )

# sampleData <- prepare_corrmat_data(
#   subset = 1:p,
#   healthy_index_name = 'CONTROLS',

#   link = "main_work/Data/Amnesia_all_AAL.mat",
#   corr_matrix_name = 'corrmats',
#   sick_index_name = 'TGA'

#   link = "main_work/Data/ADNI_data_AD_CN.mat",
#   corr_matrix_name = 'all.corrmats',
#   sick_index_name = 'AD'

#   link = "main_work/Data/NMDA_all_data_AAL90.mat",
#   corr_matrix_name = 'group.all',
#   sick_index_name = 'NMDA'
# )

test_corr_mat(sampleData)

Pelet_Cov <- with(sampleData$samples, estimate_alpha(healthy_dt = healthy, sick_dt = sick))

gee_var <- with(sampleData$samples, compute_gee_variance(
  healthy_dt = healthy, sick_dt = sick, cov_obj = Pelet_Cov
  ))

# steps <- transpose(Pelet_Cov$steps)
# steps$theta <- t(do.call(cbind, steps$theta))
# steps$alpha <- t(do.call(cbind, steps$alpha))
# steps$value <- do.call(c, steps$value)
# 
# z <- (Pelet_Cov$alpha - 1)/sqrt_diag(gee_var)
# p <- 2*pnorm(abs(z), lower.tail = F)
# p.adjust(p, 'BH') < 0.2
# p < 0.05

sd(Pelet_Cov$alpha)

mean_sqrt_diag(
  with(sampleData$samples, compute_gee_variance(
    healthy_dt = healthy, sick_dt = sick, cov_obj = Pelet_Cov, correct = T
  ))
)

mean_sqrt_diag(
  with(sampleData$samples, compute_gee_variance(
    healthy_dt = healthy, sick_dt = sick, cov_obj = Pelet_Cov, correct = F
  ))
)


# todo: estimated n?
# 
# gc()
# 
# # Pelet_Cov_jacknife <- estimateAlpha_jacknife(
# #   healthy.data = sampleData$samples$healthy, sick.data = sampleData$samples$sick,
# #   dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
# #   T_thresh = Tlength, updateU = 1, progress = T, linkFun = linkFun, jack_healthy = TRUE)
# # 
# # alpha_jk_estimate <- colMeans(Pelet_Cov_jacknife$alpha)
# # alpha_jk_variance <- var(Pelet_Cov_jacknife$alpha)*nrow(Pelet_Cov_jacknife$alpha - 1)
# 
# # save.image('main_work/data/enviroments/full_run_jacknife_simulated_data_multiplicative_link.RData')
# 
