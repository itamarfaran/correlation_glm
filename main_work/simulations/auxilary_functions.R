source('main_work/code/01_general_functions.R')
source('main_work/code/02_simulation_functions.R')
source('main_work/code/03_estimation_functions.R')
source('main_work/code/04_inference_functions.R')
ipak('gridExtra', 'latex2exp')

theme_user <- theme_bw
plot_files_path <- 'main_work/simulations/'
dpi <- 100

custom_ggsave <- function(filename, plot, width = 1, height = 1, ...){
  ggsave(filename = filename, plot = plot, path = plot_files_path,
         width = width*10, height = height*10, units = 'cm', dpi = dpi, ...)
}
  
create_estimates <- function(n_sim, n, p, percent_alpha, range_alpha, ARMA = 0, verbose = FALSE){
  # no arma, null case
  case = if(percent_alpha == 0) 'No Effect' else 'Effect'
  autocorrelated = if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if(ARMA == 0) ARMA <- NULL
  
  samples <- create_samples(n_sim = n_sim, n_h = n/2, n_s = n/2, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA)
  results <- estimate_alpha(samples$samples$healthy, samples$samples$sick, verbose=verbose)
  
  theta_dt <- data.table(type = 'Theta', Estimate = results$theta, Parameter = triangle2vector(samples$real_theta))
  alpha_dt <- data.table(type = 'Alpha', Estimate = as.vector(results$alpha), Parameter = as.vector(samples$alpha))
  
  out <- rbind(alpha_dt, theta_dt)
  out[,`:=`(case = case, autocorrelated = autocorrelated, n = n, p = p)]
  out[type == 'Alpha', Value := ifelse(Parameter == 1, 'Null', 'Non-Null')]
  return(out)
}

create_variance_estimates <- function(n_sim, n, p, p_s, percent_alpha, range_alpha, ARMA = 0){
  case = if(percent_alpha == 0) 'No Effect' else "Effect"
  autocorrelated = if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(p_s*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA)
  results <- lapply(
    1:n_sim, function(i) estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      verbose = FALSE)
  )
  
  emp_cov <- var(t(do.call(cbind, transpose(results)$alpha)))
  
  gee_vars <- calculate_mean_matrix(sapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick
  ), simplify = 'array'))
  
  out <- data.table(n = n, p = p, alpha = as.vector(samples$alpha), emp = sqrt_diag(emp_cov), est = sqrt_diag(gee_vars))
  out[,`:=`(autocorrelated = autocorrelated, p_s = p_s, case = case)]
  return(out)
}

