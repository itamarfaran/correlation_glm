source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 32
n_h <- 12
n_s <- 19
T_thresh <- Tlength <- 115
n_sim <- 70

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, check_invertability_arma)

jacknife_function <- function(i){
  estimates <- estimate_alpha_jacknife(
    healthy_dt = sample_data$samples[[i]]$healthy,
    sick_dt = sample_data$samples[[i]]$sick,
    ncores = 1, verbose = F)
  
  inference <- infer_jacknife(estimates)
  estimates[['alpha_estimate']] <- inference$estimate
  estimates[['alpha_variance']] <- inference$variance
  
  return(estimates)
}


sample_data <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = Tlength,
                              percent_alpha = 0.4, range_alpha = c(0.6, 0.8), ncores = ncores)
bootstrap_results <- pbmclapply(seq_len(n_sim), jacknife_function, mc.cores = ncores)
bootstrap_results_t <- purrr::transpose(bootstrap_results)

estimates <- do.call(rbind, bootstrap_results_t$alpha_estimate)
variances <- simplify2array(bootstrap_results_t$alpha_variance)
variances_dt <- data.table(
  empiric = triangle2vector(var(estimates)),
  theoretical = triangle2vector(calculate_mean_matrix(variances))
)
variances_diag_dt <- data.table(
  empiric = diag(var(estimates)),
  theoretical = diag(calculate_mean_matrix(variances))
)


z_values <- (estimates - linkFun$NULL_VAL)/do.call(rbind, lapply(bootstrap_results_t$alpha_variance, sqrt_diag))
p_values <- 2*pnorm(abs(z_values), lower.tail = F)
hist(z_values[,sample_data$alpha == 1])
hist(p_values[,sample_data$alpha == 1])

qplot(variances_dt$theoretical, variances_dt$empiric) + geom_smooth(method = 'lm')
qplot(variances_diag_dt$theoretical, variances_diag_dt$empiric) + geom_smooth(method = 'lm')

summary(lm(empiric ~ theoretical, variances_dt))
save.image('dev/simulations/check_jacknife.RData')
