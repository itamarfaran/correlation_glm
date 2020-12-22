source('lib/utils/general_functions.R')
source('lib/utils/preprocess_functions.R')
source('lib/utils/simulation_functions.R')
source('lib/model/estimation_functions.R')
source('lib/model/inference_functions.R')


##### simulation pars #####

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


##### data analysis configs pars #####

data_configs <- list(
  amnesia_aal = list(
    link = 'lib/data/Amnesia_all_AAL.mat',
    corr_matrix_name = 'corrmats',
    healthy_index_name = 'CONTROLS',
    sick_index_name = 'TGA'
  ),
  amnesia_aicha = list(
    link = 'lib/data/Amnesia_all_AICHA.mat',
    corr_matrix_name = 'corrmats',
    healthy_index_name = 'CONTROLS',
    sick_index_name = 'TGA'
  ),
  adni_aal = list(
    link = 'lib/data/ADNI_data_AD_CN.mat',
    corr_matrix_name = 'all.corrmats',
    healthy_index_name = 'CONTROLS',
    sick_index_name = 'AD'
  ),
  nmda_aal = list(
    link = 'lib/data/NMDA_all_data_AAL90.mat',
    corr_matrix_name = 'group.all',
    healthy_index_name = 'CONTROLS',
    sick_index_name = 'NMDA'
  )
)


##### data analysis pars #####
data_conf_name <- 'amnesia_aal'
data_conf <- if(data_conf_name == 'simulation') NULL else data_configs[[data_conf_name]]
data_subset <- NULL

if(is.null(data_conf)){
  sample_data <- create_samples(
    n_h = n_h, n_s = n_s, p = p, Tlength = Tlength, dim_alpha = 1,
    percent_alpha = 0.3, range_alpha = c(0.9, 0.9),
    ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
    ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
    ncores = ncores
  )
} else {
  sample_data <- prepare_corrmat_data(
    link = data_conf$link,
    corr_matrix_name = data_conf$corr_matrix_name,
    healthy_index_name = data_conf$healthy_index_name,
    sick_index_name = data_conf$sick_index_name,
    subset = if(is.null(data_subset)) NULL else 1:data_subset
  )
}

##### analysis #####

test_corr_mat(sample_data)

results <- estimate_alpha(
  healthy_dt = sample_data$samples$healthy,
  sick_dt = sample_data$samples$sick,
  early_stop = FALSE,
  matrix_reg_config = list(method = 'constant', const = 0.1)
  )

gee_var <- compute_gee_variance(
  healthy_dt = sample_data$samples$healthy,
  sick_dt = sample_data$samples$sick,
  cov_obj = results
  )

# results_jacknife <- estimate_alpha_jacknife(
#   healthy_dt = sample_data$samples$healthy, sick_dt = sample_data$samples$sick,
#   linkFun = linkFun, jack_healthy = TRUE, return_gee = F, ncores = ncores)
# jacknife_inference <- infer_jacknife(results_jacknife)

save.image(paste0('lib/data/envs/full_run_jacknife_', data_conf_name, '_multiplicative_link.RData'))

