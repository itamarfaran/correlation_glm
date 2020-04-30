source("main_work/code/01_general_functions.R")
source("main_work/code/02_simulation_functions.R")
source("main_work/code/03_estimation_functions.R")
source("main_work/code/04_inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
B <- 200

bootstrap_gee <- function(p, n, percent_alpha = 0.3, range_alpha = c(1, 1),
                          sick_obs_percentage = 0.5, t_length = 115,
                          ar = 0, nboot = 100, sig_level = 0.05,
                          sim_num = NULL, ncores = 1, verbose = FALSE){
  lapply_ <- if(verbose) pbmcapply::pbmclapply else parallel::mclapply

  samples <- create_samples(
    n_sim = nboot,
    n_h = round(n*(1 - sick_obs_percentage)),
    n_s = round(n*sick_obs_percentage),
    p = p,
    Tlength = t_length,
    dim_alpha = 1,
    percent_alpha = percent_alpha, range_alpha = range_alpha,
    ARsick = ar, MAsick = NULL,
    ARhealth = ar, MAhealth = NULL,
    linkFun = linkFun,
    ncores = ncores
  )
  
  estimates <- lapply_(1:nboot, function(i){
    sample_ <- samples$samples[[i]]
    
    covobj <- estimate_alpha(
      healthy_dt = sample_$healthy,
      sick_dt = sample_$sick,
      dim_alpha = 1, verbose = F, linkFun = linkFun)
    
    covobj$gee_var_new <- compute_gee_variance(
      covobj, 
      healthy_dt = sample_$healthy,
      sick_dt = sample_$sick
    )
    covobj$steps <- covobj$Log_Optim <- NULL
    return(covobj)
  }, mc.cores = ncores)
  
  alpha_mat <- do.call(rbind, lapply(transpose(estimates)$alpha, as.vector))
  alpha_sd_mat <- do.call(rbind, lapply(transpose(estimates)$gee_var_new, sqrt_diag))
  alpha_real <- as.vector(samples$alpha)

  out <- data.table(
    sim_num = sim_num,
    p = p,
    n = n,
    percent_alpha = percent_alpha,
    min_alpha = min(range_alpha),
    sick_obs_percentage = sick_obs_percentage,
    Tlength = t_length,
    ar = ar,
    confidence = 1 - sig_level,
    alpha = alpha_real,
    alpha_est_mean = colMeans(alpha_mat),
    # alpha_est_emp_sd = apply(alpha_mat, 2, sd_known_mu, mu = alpha_real),
    alpha_est_emp_sd = apply(alpha_mat, 2, sd),
    # alpha_est_lower = apply(alpha_mat, 2, quantile, probs = sig_level/2),
    # alpha_est_upper = apply(alpha_mat, 2, quantile, probs = 1 - sig_level/2),
    gee_sd_mean = colMeans(alpha_sd_mat),
    gee_sd_emp_sd = apply(alpha_sd_mat, 2, sd)
    # gee_sd_lower = apply(alpha_sd_mat, 2, quantile, probs = sig_level/2),
    # gee_sd_upper = apply(alpha_sd_mat, 2, quantile, probs = 1 - sig_level/2)
  )
  return(out)
}

combinations2boot <- data.table(
  p = 10*round((40*rbeta(B, 1, 1.8) + 20)/10),
  n = 10*round(runif(B, 80, 120)/10),
  sick_obs_percentage = 10*round((0.6*rbeta(B, 1.5, 1.5) + 0.2)/10, 2),
  t_length = 10*round((40*rbeta(B, 1.5, 1.5) + 80)/10),
  percent_alpha = sample(c(0, 0.1, 0.2, 0.3), B, TRUE),
  min_alpha = sample(c(0.7, 0.8, 0.9), B, TRUE),
  ar = round(c(
    runif(B/4, -0.7, -0.3),
    rep(0, B/2),
    runif(B/4, 0.3, 0.7)
  ), 1)
)
combinations2boot[percent_alpha == 0, min_alpha := 1]

results <- do.call(rbind, pbmclapply(1:B, function(b){
  with(combinations2boot[b], bootstrap_gee(
    p = p,
    n = n,
    sick_obs_percentage = sick_obs_percentage,
    t_length = t_length,
    percent_alpha = percent_alpha,
    range_alpha = c(min_alpha, 1),
    ar = ar,
    nboot = ncores,
    ncores = ncores
    ))
  }, mc.cores = 1))

fwrite(results, 'main_work/code/gee-bootstrap-app/gee_data_boot.csv')
save.image(paste0('main_work/data/enviroments/gee_new_bootstrap', format(Sys.time(), '%Y%m%d_%H%M'), '.RData'))
