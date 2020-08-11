##### source + functions #####

source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")


#  todo: check in high and low sample size in comparison to the REAL alpha
#  todo: Show path with regularization
#  todo: estimate OLS, Regularized, and Full on Random/Non Random Effect (total 6)
#  check best estimators MSE and Bias; Simulate with real TGA Theta

get_par_path <- function(model, par = 'alpha', model_name = NULL){
  transpose(model$steps)[[par]] %>% 
    lapply(as.vector) %>% do.call(cbind, .) %>%
    data.table() -> dt
  
  colnames(dt) <- as.character(seq_len(ncol(dt)))
  dt[,index := 1:.N]
  dt <- melt(dt, id.vars = 'index', variable.name = 'step', value.name = 'parameter')
  dt[,step := as.numeric(step)]
  dt[,index := factor(index)]
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

ggsave_batch <- function(l_plots, names, dir = 'temp'){
  if(!dir.exists(dir)) dir.create(dir)
  for(i in seq_along(l_plots)) ggsave(paste0(dir, '/', names[i], '.png'), l_plots[[i]])
}

##### definitions #####
p <- 60
linkFun <- linkFunctions$multiplicative_identity
sample_size_plt <- 20

##### estimate on tga data #####
tga_data <- prepare_corrmat_data(
  subset = 1:p,
  healthy_index_name = 'CONTROLS',
  link = "lib/data/Amnesia_all_AAL.mat",
  corr_matrix_name = 'corrmats',
  sick_index_name = 'TGA'
)

tga_iid_model <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$sick),
  dim_alpha = 1, linkFun = linkFun,
  cov_method = 'identity'
)

tga_cov_model <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$sick),
  alpha0 = tga_iid_model$alpha, theta0 = tga_iid_model$theta,
  linkFun = tga_iid_model$linkFun,
  # matrix_reg_config = list(do_reg = T, method = 'increase_diag', const = 0.75),
  cov_method = 'corrmat'
)

tga_cov_model_reg <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$sick),
  alpha0 = tga_iid_model$alpha, theta0 = tga_iid_model$theta,
  linkFun = tga_iid_model$linkFun,
  matrix_reg_config = list(do_reg = T, method = 'increase_diag', const = 0.25),
  cov_method = 'corrmat'
)

tga_cov_model_start_null <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(tga_data$samples$sick),
  linkFun = tga_iid_model$linkFun,
  # matrix_reg_config = list(do_reg = T, method = 'increase_diag', const = 0.25),
  cov_method = 'corrmat'
)

##### plot tga data #####

tga_alpha_steps <- get_full_path(
  get_par_path(tga_iid_model, 'alpha', 'iid'),
  get_par_path(tga_cov_model, 'alpha', 'cov')
)
tga_theta_steps <- get_full_path(
  get_par_path(tga_iid_model, 'theta', 'iid'),
  get_par_path(tga_cov_model, 'theta', 'cov')
)

set.seed(847)
plot_list <- list(
  create_plot(tga_alpha_steps, hline = 1, size = sample_size_plt),
  create_plot(tga_theta_steps, hline = 0, size = sample_size_plt),
  qplot(x = tga_iid_model$alpha, y = tga_cov_model$alpha) +
    geom_abline(slope = 1, intercept = 0),
  qplot(x = tga_cov_model$alpha, y = tga_cov_model_reg$alpha) +
    geom_abline(slope = 1, intercept = 0),
  qplot(x = tga_cov_model$alpha, y = tga_cov_model_start_null$alpha) +
    geom_abline(slope = 1, intercept = 0),
  qplot(x = tga_iid_model$alpha, y = tga_cov_model_reg$alpha) +
    geom_abline(slope = 1, intercept = 0)
)

tga_mcd_results <- minimize_cov_distance(tga_data$samples$sick, tga_iid_model)

tga_alpha_steps_reg <- get_full_path(
  get_par_path(tga_iid_model, 'alpha', 'iid'),
  get_par_path(tga_cov_model_reg, 'alpha', 'cov')
)

set.seed(847)
plot_list[[length(plot_list) + 1]] <- create_plot(tga_alpha_steps_reg, hline = 1, size = sample_size_plt)
ggsave_batch(plot_list, c('tga_alpha', 'tga_theta', 'tga_iid_vs_cov', 'tga_cov_reg',
                          'tga_cov_start_points', 'tga_alpha_reg'))

##### estimate on simulated data, high sample size #####
ARMAdetails <- list(
  ARsick = c(0.5, 0.1), MAsick = c(0.5, 0.1),
  ARhealth = c(0.4, 0.2), MAhealth = c(0.4, 0.2)
)
percent_alpha <- 0.2
range_alpha <- c(0.85, 1.05)


n_h <- 180
n_s <- 180
T_thresh <- Tlength <- 120

high_sample_data <- create_samples(
  n_h = n_h, n_s = n_s, p = p, Tlength = Tlength, dim_alpha = 1,
  percent_alpha = percent_alpha, range_alpha = range_alpha,
  ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
  ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
  ncores = ncores
)

high_iid_model <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(high_sample_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(high_sample_data$samples$sick),
  dim_alpha = 1, linkFun = linkFun,
  cov_method = 'identity'
)

high_cov_model <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(high_sample_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(high_sample_data$samples$sick),
  alpha0 = high_iid_model$alpha, theta0 = high_iid_model$theta,
  linkFun = high_iid_model$linkFun,
  # matrix_reg_config = list(do_reg = T, method = 'increase_diag', const = 0.75),
  cov_method = 'corrmat'
)

##### plot simulated data #####

high_alpha_steps <- get_full_path(
  get_par_path(high_iid_model, 'alpha', 'iid'),
  get_par_path(high_cov_model, 'alpha', 'cov')
)
high_theta_steps <- get_full_path(
  get_par_path(high_iid_model, 'theta', 'iid'),
  get_par_path(high_cov_model, 'theta', 'cov')
)

ggsave_batch(list(
  create_plot(high_alpha_steps, hline = 1, size = sample_size_plt),
  create_plot(high_theta_steps, hline = 0, size = sample_size_plt),
  qplot(x = high_iid_model$alpha, y = high_cov_model$alpha) +
    geom_abline(slope = 1, intercept = 0)
  ), c('high_samp_sim_alpha', 'high_samp_sim_theta', 'high_samp_sim_iid_vs_cov'))


high_mcd_results <- minimize_cov_distance(high_sample_data$samples$sick, high_iid_model)


##### estimate on simulated data, low sample size #####
n_h <- 20
n_s <- 20
T_thresh <- Tlength <- 70

low_sample_data <- create_samples(
  n_h = n_h, n_s = n_s, p = p, Tlength = Tlength, dim_alpha = 1,
  percent_alpha = percent_alpha, range_alpha = range_alpha,
  ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
  ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
  ncores = ncores
)

low_iid_model <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(low_sample_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(low_sample_data$samples$sick),
  dim_alpha = 1, linkFun = linkFun,
  cov_method = 'identity'
)

low_cov_model <- estimate_loop(
  healthy_dt = convert_corr_array_to_data_matrix_test(low_sample_data$samples$healthy),
  sick_dt = convert_corr_array_to_data_matrix_test(low_sample_data$samples$sick),
  alpha0 = low_iid_model$alpha, theta0 = low_iid_model$theta,
  linkFun = low_iid_model$linkFun,
  # matrix_reg_config = list(do_reg = T, method = 'increase_diag', const = 0.75),
  cov_method = 'corrmat'
)

##### plot simulated data #####

low_alpha_steps <- get_full_path(
  get_par_path(low_iid_model, 'alpha', 'iid'),
  get_par_path(low_cov_model, 'alpha', 'cov')
)


ggsave_batch(list(
  create_plot(low_alpha_steps, hline = 1, size = sample_size_plt),
  qplot(x = low_iid_model$alpha, y = low_cov_model$alpha) +
    geom_abline(slope = 1, intercept = 0)
), c('low_samp_sim_alpha', 'low_samp_sim_iid_vs_cov'))


low_mcd_results <- minimize_cov_distance(low_sample_data$samples$sick, low_iid_model)

##### end #####
warnings()

save.image('temp/check_main_code.RData')

message('tga mcd')
print(tga_mcd_results)

message('high sample mcd')
print(high_mcd_results)

message('low sample mcd')
print(low_mcd_results)
