source("tex/simulations/aux_.R")

linkFun <- LinkFunctions$multiplicative_identity


conf <- list(
  TGA = list(
    link = "lib/data/Amnesia_all_AAL.mat",
    corr_matrix_name = 'corrmats',
    sick_index_name = 'TGA'
  )
)
desease_data <- 'TGA'
file <- paste0('tex/tga_analysis/analysis_', tolower(desease_data), '_', linkFun$NAME, 'nobiascorr.RData')

if (file.exists(file)){
  load(file)
} else {
  sample_data <- prepare_corrmat_data(
    healthy_index_name = 'CONTROLS',
    link = conf[[desease_data]]$link,
    corr_matrix_name = conf[[desease_data]]$corr_matrix_name,
    sick_index_name = conf[[desease_data]]$sick_index_name
  )
  
  sapply(sample_data$samples, test_corr_mat)
  
  results <- estimate_model(
    control_arr = sample_data$samples$healthy,
    diagnosed_arr = sample_data$samples$sick,
    LinkFunc = linkFun,
    bias_correction = FALSE)
  
  gee_var <- with(sample_data$samples, compute_gee_variance(
    control_arr = healthy, diagnosed_arr = sick, mod = results, est_mu = T
  ))
  
  save(sample_data, results, gee_var, file = file)
}

out <- data.table(
  desease_data = desease_data,
  estimate = as.vector(results$alpha),
  sd = sqrt_diag(gee_var)
)
out[,z_value := (estimate - linkFun$null_value)/sd]
out[,p_value := 2*pnorm(abs(z_value), lower.tail = F)]
out[,p_adjusted := p.adjust(p_value, 'BH')]
out[,index := 1:.N]
print(out[p_adjusted < .05])
