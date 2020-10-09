source("lib/utils/general_functions.R")
source("lib/utils/preprocess_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")
ipak('gridExtra', 'latex2exp')

theme_user <- theme_bw
plot_files_path <- 'tex/simulations/' 
dpi <- 100

reverselog_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(
    paste0("reverselog-", format(base)), trans, inv,
    scales::log_breaks(base = base),
    domain = c(1e-100, Inf))
}

square_plot <- function(plt){
  r <- max(abs(layer_scales(plt)$x$range$range))
  s <- max(abs(layer_scales(plt)$y$range$range))
  t <- round(max(r, s), 1)
  plt <- plt + coord_equal(xlim = c(0, t), ylim = c(0, t))
  return(plt)
}

custom_ggsave <- function(filename, plot, width = 1, height = 1, ...){
  ggsave(filename = filename, plot = plot, path = plot_files_path,
         width = width*10, height = height*10, units = 'cm', dpi = dpi, ...)
}

ggsave_batch <- function(l_plots, names, dir = 'temp'){
  if(!dir.exists(dir)) dir.create(dir)
  for(i in seq_along(l_plots)) ggsave(paste0(dir, '/', names[i], '.png'), l_plots[[i]])
}

get_par_path <- function(model, par = 'alpha', model_name = NULL, real_par = NULL){
  transpose(model$steps)[[par]] %>%
    lapply(as.vector) %>% do.call(cbind, .) %>%
    data.table() -> dt

  colnames(dt) <- as.character(seq_len(ncol(dt)))
  dt[,index := 1:.N]
  dt <- melt(dt, id.vars = 'index', variable.name = 'step', value.name = 'parameter')
  dt[,step := as.numeric(step)]
  dt[,index := factor(index)]
  if(!is.null(real_par)) dt[,real_par := rep(real_par, times = max(step))]
  if(!is.null(model_name)) dt[,model := model_name]
  return(dt)
}

get_full_path <- function(...){
  dt_list <- list(...)
  for(i in 2:length(dt_list)) dt_list[[i]][,step := step + max(dt_list[[i - 1]]$step) - 1]
  dt_full <- do.call(rbind, dt_list)
  return(dt_full)
}

create_plot <- function(dt, hline = NULL, size = 20){
  if(size < 1 | size > nrow(dt)) size <- nrow(dt)
  samp <- dt[,sort(sample(unique(index), size))]
  p <- ggplot(dt[index %in% samp], aes(x = step, y = parameter, col = index))
  if(!is.null(hline)) p <- p + geom_hline(yintercept = hline)
  p <- p + geom_vline(xintercept = dt[model == 'iid', max(step)],
                      color = 'lightgrey', linetype = 2, size = 1) +
    geom_line() + geom_point(aes(shape = model)) +
    theme_bw() +
    theme(legend.position = 'none')
  if('real_par' %in% colnames(dt)) p <- p + geom_hline(aes(yintercept = real_par, col = index), linetype = 2)

  return(p)
}

get_emp_cov <- function(dt, model){
  expected_dt <- rep(1, nrow(dt)) %o% triangle2vector(with(model, linkFun$FUN(theta, alpha, 1)))
  residuals_ <- dt - expected_dt
  emp_cov <- t(residuals_) %*% residuals_

  return(emp_cov)
}

minimize_cov_distance <- function(data, model, reg_method = 'increase_diag'){
  reg_method <- match.arg(reg_method, c('avg_diag', 'increase_diag'), F)
  dt <- convert_corr_array_to_data_matrix_test(data)
  emp_cov <- get_emp_cov(dt, model)
  g12 <- corrmat_covariance_from_dt(dt)

  optim_fun <- function(par, reg_method) sum(triangle2vector(
    regularize_matrix(g12, reg_method, par, F) - emp_cov
  )^2)

  out <- optimise(optim_fun, 0:1, reg_method = reg_method)
  return(out)
}

create_estimates <- function(n_sim, n, p, percent_alpha, range_alpha, ARMA = 0, verbose = FALSE){
  # no arma, null case
  case <- if(percent_alpha == 0) 'No Effect' else 'Effect'
  autocorrelated <- if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
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

create_variance_estimates <- function(n_sim, n, p, p_s, percent_alpha, range_alpha, ARMA = 0, ncores = 1){
  case <- if(percent_alpha == 0) 'No Effect' else "Effect"
  autocorrelated <- if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(p_s*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA, ncores = ncores)
  results <- pbmclapply(
    1:n_sim, function(i) estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      verbose = FALSE), mc.cores = ncores
  )
  
  emp_cov <- var(t(do.call(cbind, transpose(results)$alpha)))
  
  gee_vars <- calculate_mean_matrix(simplify2array(pbmclapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick
  ), mc.cores = ncores)))
  
  out <- data.table(n = n, p = p, alpha = as.vector(samples$alpha), emp = sqrt_diag(emp_cov), est = sqrt_diag(gee_vars))
  out[,`:=`(autocorrelated = autocorrelated, p_s = p_s, case = case)]
  return(out)
}

create_sample_estimates <- function(n_sim, n, p, p_s, percent_alpha, range_alpha, ARMA = 0, ncores = 1, sim = NULL){
  case <- if(percent_alpha == 0) 'No Effect' else "Effect"
  autocorrelated <- if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(p_s*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA,
                            ncores = ncores)
  results <- pbmclapply(
    1:n_sim, function(i) estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      verbose = FALSE), mc.cores = ncores
  )
  
  gee_vars <- pbmclapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick
  ), mc.cores = ncores)
  
  out <- data.table(
    sim_num = rep(1:n_sim, each = p),
    voxel = rep(1:p, times = n_sim),
    estimated_alpha = do.call(c, transpose(results)$alpha),
    real_alpha = rep(samples$alpha, times = n_sim),
    sd = do.call(c, lapply(gee_vars, sqrt_diag))
  )
  out[,`:=`(autocorrelated = autocorrelated, p_s = p_s, case = case)]
  if(!is.null(sim)) out[,sim := sim]
  return(out)
}

