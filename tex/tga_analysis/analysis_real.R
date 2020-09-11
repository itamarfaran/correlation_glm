source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity


conf <- list(
  TGA = list(
    link = "lib/data/Amnesia_all_AAL.mat",
    corr_matrix_name = 'corrmats',
    sick_index_name = 'TGA'
  ),
  AD = list(
    link = "lib/data/ADNI_data_AD_CN.mat",
    corr_matrix_name = 'all.corrmats',
    sick_index_name = 'AD'
  ),
  NMDA = list(
    link = "lib/data/NMDA_all_data_AAL90.mat",
    corr_matrix_name = 'group.all',
    sick_index_name = 'NMDA'
  )
)
desease_data <- 'TGA'
file <- paste0('tex/tga_analysis/analysis_', tolower(desease_data), '_', linkFun$NAME, '.RData')

if (file.exists(file)){
  load(file)
} else {
  sample_data <- prepare_corrmat_data(
    healthy_index_name = 'CONTROLS',
    link = conf[[desease_data]]$link,
    corr_matrix_name = conf[[desease_data]]$corr_matrix_name,
    sick_index_name = conf[[desease_data]]$sick_index_name
  )
  
  test_corr_mat(sample_data)
  
  results <- with(sample_data$samples, estimate_alpha(healthy_dt = healthy, sick_dt = sick, linkFun = linkFun))
  
  gee_var <- with(sample_data$samples, compute_gee_variance(
    healthy_dt = healthy, sick_dt = sick, cov_obj = results, est_mu = F
  ))
  
  save(sample_data, results, gee_var, file = file)
}

out <- data.table(
  desease_data = desease_data,
  estimate = as.vector(results$alpha),
  sd = sqrt_diag(gee_var)
)
out[,z_value := (estimate - linkFun$NULL_VAL)/sd]
out[,p_value := 2*pnorm(abs(z_value), lower.tail = F)]
out[,p_adjusted := p.adjust(p_value, 'BH')]
out[,index := 1:.N]
print(out[p_adjusted < .05])
