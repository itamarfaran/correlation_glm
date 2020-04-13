source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 25
nH <- 30
nS <- 30
T_thresh <- Tlength <- 115

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, checkInv)

sampleData <- createSamples(
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

Pelet_Cov <- estimateAlpha(
  healthy.data = sampleData$samples$healthy, sick.data = sampleData$samples$sick,
  dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
  T_thresh = T_thresh, updateU = 1, progress = T, linkFun = linkFun)

gee_var <- compute_gee_variance(
  CovObj = Pelet_Cov, sampledata = sampleData$samples,
  est_mu = TRUE, correct = FALSE
)

gc()

# Pelet_Cov_jacknife <- estimateAlpha_jacknife(
#   healthy.data = sampleData$samples$healthy, sick.data = sampleData$samples$sick,
#   dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
#   T_thresh = Tlength, updateU = 1, progress = T, linkFun = linkFun, jack_healthy = TRUE)
# 
# alpha_jk_estimate <- colMeans(Pelet_Cov_jacknife$alpha)
# alpha_jk_variance <- var(Pelet_Cov_jacknife$alpha)*nrow(Pelet_Cov_jacknife$alpha - 1)

# save.image('main_work/data/enviroments/full_run_jacknife_simulated_data_multiplicative_link.RData')

