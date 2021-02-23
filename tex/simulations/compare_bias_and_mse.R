source("tex/simulations/aux_.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 32
n_h <- 12
n_s <- 19
T_thresh <- Tlength <- 115
n_sim <- 50
early_stop <- TRUE

ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
sapply(ARMAdetails, is_invertable_arma)

patient_data <- prepare_corrmat_data(
  subset = 1:p,
  healthy_index_name = 'CONTROLS',
  link = "lib/data/Amnesia_all_AAL.mat",
  corr_matrix_name = 'corrmats',
  sick_index_name = 'TGA'
)
real_theta <- calculate_mean_matrix(patient_data$samples$healthy)

sample_data <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = Tlength,
                              percent_alpha = 0.4, range_alpha = c(0.6, 0.8), real_theta = real_theta, ncores = ncores)

estimate_fun <- function(i){
  healthy_dt <- convert_corr_array_to_data_matrix(sample_data$samples[[i]]$healthy)
  sick_dt <- convert_corr_array_to_data_matrix(sample_data$samples[[i]]$sick)
  weight_matrix <- corrmat_covariance_from_datamatrix(sick_dt)
  
  ols <- corrfuncs:::optimiser(healthy_dt = healthy_dt, sick_dt = sick_dt, weight_matrix = NULL, early_stop = early_stop, verbose = F)
  wls_warm_start <- inner_optim_loop(healthy_dt = healthy_dt, sick_dt = sick_dt, weight_matrix = weight_matrix, early_stop = early_stop,
                                     alpha0 = ols$alpha, theta0 = ols$theta, verbose = F)
  wls <- corrfuncs:::optimiser(healthy_dt = healthy_dt, sick_dt = sick_dt, weight_matrix = weight_matrix, early_stop = early_stop, verbose = F)
  wls_reg <- corrfuncs:::optimiser(healthy_dt = healthy_dt, sick_dt = sick_dt, weight_matrix = weight_matrix, early_stop = early_stop,
                              matrix_reg_config = list(method = 'increase_diag', const = 0.25), verbose = F)
  
  estimates <- data.table(
    sim_num = i,
    ind = 1:p,
    ols = as.vector(ols$alpha),
    wls = as.vector(wls$alpha),
    wls_reg = as.vector(wls_reg$alpha),
    wls_warm_start = as.vector(wls_warm_start$alpha)
  )
  
  to_bind <-list(
    ols$convergence,
    wls$convergence,
    wls_warm_start$convergence,
    wls_reg$convergence
  )
  max_len <- max(sapply(to_bind, length))
  for (j in seq_along(to_bind)) to_bind[[j]] <- c(to_bind[[j]], rep(NA, max_len - length(to_bind[[j]]) ))
  to_bind <- as.data.table(to_bind)
  to_bind <- to_bind[2:.N]
  to_bind[,sim_num := i]
  colnames(to_bind) <- c('ols', 'wls', 'wls_reg', 'wls_warm_start', 'sim_num')
  setcolorder(to_bind, 'sim_num')
  
  return(list(estimates = estimates, convergence = to_bind))
}

results <- pbmclapply(seq_len(n_sim), estimate_fun, mc.cores = ncores)

convergence <- rbindlist(purrr::transpose(results)$convergence)

estimates <- rbindlist(purrr::transpose(results)$estimates)
estimates[,alpha := rep(sample_data$alpha, n_sim)]
estimates <- melt(estimates, id.vars = c('sim_num', 'ind', 'alpha'), variable.name = 'method', value.name = 'estimate')

bias_dt <- estimates[,.(bias = mean(estimate - alpha), rmse = sqrt(mean((estimate - alpha)^2))), by = .(method, ind)]
bias_dt <- melt(bias_dt, id.vars = c('method', 'ind'), variable.name = 'type', value.name = 'value')

bias_plot <- ggplot(bias_dt[type == 'bias'], aes(x = value, fill = method)) +
  geom_histogram(color = 'black') + facet_grid(method ~ .) + 
  geom_vline(xintercept = 0, size = 2) + 
  theme_user() + theme(legend.position = 'none')

rmse_plot <- ggplot(bias_dt[type == 'rmse'], aes(x = value, fill = method)) +
  geom_histogram(color = 'black') + facet_grid(method ~ .) + 
  theme_user() + theme(legend.position = 'none')

means <- bias_dt[,.(means = mean(value), medians = median(value)), by=.(method, type)][order(means)]

save.image('tex/simulations/compare_bias_and_mse.RData')

