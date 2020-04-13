##### config #####

B = 100
n_sample = 50
n_repeat = 5
ar = c(0.2, 0.5)
percent_sick = 0.5
theta = runif(n_repeat)
alpha = rep(1, n_repeat) # runif(n_repeat)

ncores <<- max(parallel::detectCores() - 2, 1)
do_save <- FALSE

##### functions & ipak #####

ipak <- function(..., only_install = FALSE){
  pkg <- unlist(list(...))
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  if(!only_install) sapply(pkg, require, character.only = TRUE)
}

ipak('data.table', 'dplyr', 'purrr', 'stringr', 'geepack', 'numDeriv', 'ggplot2', 'ParallelLogger')

sqrt_diag <- function(x) sqrt(diag(x))
simplify2array_ <- function(x, higher = TRUE){
  out <- simplify2array(x = x, higher = higher)
  if(is.matrix(out)) out <- t(out)
  return(out)
}
extract_number <- function(x) as.numeric(str_extract(x, "(\\d)+"))
create_ar_matrix <- function(p, ar){
  out <- matrix(0, p, p)
  out <- ar^abs(row(out) - col(out))
  return(out)
}
create_ar_vector <- function(n, ar = 0, sd = 1){
  x <- numeric(n)
  x[1] <- rnorm(1, 0, sd)
  for(i in 2:n) x[i] <- ar*x[i - 1] + rnorm(1, 0, sd)
  return(x)
}
create_gee_data_simple <- function(
  n_sample, n_repeat, ar,
  percent_sick = 0.5, 
  theta = 1, alpha = 1,
  sd = .1, sim_num = 1, seed_x = NULL){
  
  if(length(theta) == 1) theta <- rep(theta, n_repeat)
  if(length(alpha) == 1) alpha <- rep(alpha, n_repeat)
  if(length(ar) <= 1) ar <- rep(ar, 2)
  
  sick_index <- c(
    rep(0, round(n_sample*(1 - percent_sick))),
    rep(1, round(n_sample*percent_sick))
  )
  
  theta_matrix <- rep(1, n_sample) %o% theta
  alpha_matrix <- rep(1, n_sample) %o% alpha
  is_sick_matrix <- sick_index %o% rep(1, n_repeat)
  # expected_values <- theta_matrix * (alpha_matrix ^ is_sick_matrix)
  expected_values <- theta_matrix + (alpha_matrix * is_sick_matrix)
  ar_errors <- t(sapply(
      1:n_sample,
      function(i) create_ar_vector(n_repeat, ar[sick_index[i] + 1], sd = sd)
    )
  )
  
  y_lin <- expected_values + ar_errors
  colnames(y_lin) <- paste0('response_', 1:n_repeat)
  
  sample_out <- data.table(
    simulation = sim_num,
    id = 1:n_sample,
    is_sick = sick_index,
    y_lin
  )
  
  sample_out <- sample_out %>% 
    melt(
      id.vars = c('simulation', 'id', 'is_sick'),
      variable.name = 'response_index',
      value.name = 'response_value') %>% 
    mutate(response_index = extract_number(response_index)) %>%
    arrange(id, response_index) %>%
    as.data.table()
  
  summary_out <- data.table(
    simulation = sim_num,
    n_sample = n_sample,
    n_repeat = n_repeat,
    ar_healthy = ar[1],
    ar_sick = ar[2],
    sd = sd
  )
  
  parameters_out <- data.table(
    simulation = sim_num,
    voxel = 1:n_repeat,
    theta = theta, alpha = alpha
  )
  
  return(list(
    sample = sample_out,
    summary = summary_out,
    parameters = parameters_out
  ))
}


aggregate_gee_data_simple <- function(gee_data_simple){
  reverse_list <- purrr::transpose(gee_data_simple)
  
  samples <- do.call(rbind, reverse_list$sample)
  summaries <- do.call(rbind, reverse_list$summary)
  parameters <- do.call(rbind, reverse_list$parameters)
  
  return(list(
    samples = samples,
    summaries = summaries,
    parameters = parameters
  ))
}


theta_of_alpha <- function(A, healthy_dt, sick_dt){
  colMeans(rbind(
    # sick_dt / (rep(1, nrow(sick_dt)) %o% A),
    sick_dt - (rep(1, nrow(sick_dt)) %o% A),
    healthy_dt
  ))
}

split_data <- function(data, b){
  to_split <- copy(data[simulation == b])

  out <- dcast(
    to_split, 'simulation + is_sick + id ~ response_index',
    value.var = 'response_value')
  
  out <- lapply(
    split(out, by = 'is_sick'),
    function(dt) as.matrix(dt[,-c('simulation', 'is_sick', 'id')])
  )
  
  names(out) <- c('healthy', 'sick')
  return(out)
}

estimate_gee <- function(data, ar, max_iter = 250, tol = 1e-05, verbose = FALSE){
  min_function <- function(A){
    # residuals <- data$sick - (rep(1, nrow(data$sick)) %o% (theta_*A))
    residuals <- data$sick - (rep(1, nrow(data$sick)) %o% (theta_ + A))
    sigma <- create_ar_matrix(ncol(data$sick), ar[2])
    mean(diag(residuals %*% solve(sigma) %*% t(residuals)))
  }
  
  if(length(ar) <= 1) ar <- rep(ar, 2)
  
  theta_ <- colMeans(data$healthy)
  alpha_ <- optim(rep(1, ncol(data$sick)), min_function, method = 'BFGS')$par
  for(i in 1:max_iter){
    theta_0 <- theta_
    alpha_0 <- alpha_
    
    theta_ <- theta_of_alpha(A = alpha_, healthy_dt = data$healthy, sick_dt = data$sick)
    alpha_ <- optim(alpha_, min_function, method = 'BFGS')$par
    
    converg_diff_norm <- sqrt(mean((alpha_ - alpha_0)^2))
    if(verbose) message(paste0('iter: ', i, ' || diff: ', converg_diff_norm))
    if(converg_diff_norm < tol) break()
  }
  return(list(theta = theta_, alpha = alpha_, converg = 1*(i >= max_iter)))
}

compute_mu_alpha_jacobian <- function(type, alpha, healthy_dt, sick_dt){
  func <- if(type == 'sick'){
    # function(A) theta_of_alpha(A, healthy_dt = healthy_dt, sick_dt = sick_dt)*A
    function(A) theta_of_alpha(A, healthy_dt = healthy_dt, sick_dt = sick_dt) + A
  } else if(type == 'healthy') {
    function(A) theta_of_alpha(A, healthy_dt = healthy_dt, sick_dt = sick_dt) 
  }
  return(
    jacobian(func = func, x = alpha)
  )
}

compute_gee_variance <- function(estimates, data, ar, est_mu = TRUE, ncores = 1){
  
  compute_gee_raw <- function(type, list_){
    if(type == 'I0'){
      out <- t(list_$jacobian) %*% list_$solve_Sigma %*% list_$jacobian
    } else if (type == 'I1'){
      residuals <- list_$data - rep(1, nrow(list_$data)) %o% list_$expected_value
      cov_mat <- t(residuals) %*% residuals / nrow(list_$data)
      out <- t(list_$jacobian) %*% list_$solve_Sigma %*% cov_mat %*% list_$solve_Sigma %*%list_$jacobian
    }
    out <- out*nrow(list_$data)
    return(out)
  }
  
  if(length(ar) <= 1) ar <- rep(ar, 2)
  
  healthy_data <- data$healthy
  sick_data <- data$sick
  
  healthy_list <- list(
    data = healthy_data,
    jacobian = compute_mu_alpha_jacobian(
      type = 'healthy',
      alpha = estimates$alpha,
      healthy_dt = healthy_data,
      sick_dt = sick_data
    ),
    expected_value = if(est_mu) estimates$theta else colMeans(healthy_data),
    solve_Sigma = solve(create_ar_matrix(ncol(healthy_data), ar[1]))
  )
  
  sick_list <- list(
    data = sick_data,
    jacobian = compute_mu_alpha_jacobian(
      type = 'sick',
      alpha = estimates$alpha,
      healthy_dt = healthy_data,
      sick_dt = sick_data
    ),
    # expected_value = if(est_mu) estimates$theta*estimates$alpha else colMeans(sick_data),
    expected_value = if(est_mu) estimates$theta + estimates$alpha else colMeans(sick_data),
    solve_Sigma = solve(create_ar_matrix(ncol(sick_data), ar[2]))
  )
  
  I0 <- compute_gee_raw('I0', healthy_list) + compute_gee_raw('I0', sick_list)
  solve_I0 <- solve(I0)
  I1 <- compute_gee_raw('I1', healthy_list) + compute_gee_raw('I1', sick_list)
  res <- solve_I0 %*% I1 %*% solve_I0
  return(res)
}

extact_results <- function(result_list){
  result_list_trans <- transpose(result_list)
  theta_mat <- simplify2array_(result_list_trans$theta)
  alpha_mat <- simplify2array_(result_list_trans$alpha)
  convergence_vect <- simplify2array_(result_list_trans$convergence)
  vcov_arr <- simplify2array_(result_list_trans$vcov)
  sd_mat <- simplify2array_(lapply(result_list_trans$vcov, sqrt_diag))
  
  return(list(
    theta_mat = theta_mat,
    alpha_mat = alpha_mat,
    convergence_vect = convergence_vect,
    vcov_arr = vcov_arr,
    sd_mat = sd_mat
  ))
}

##### boot #####

sampling_results <- aggregate_gee_data_simple(lapply(
  1:B,
  function(b){
    create_gee_data_simple(
      n_sample = n_sample,
      n_repeat = n_repeat, 
      ar = ar,
      percent_sick = percent_sick,
      theta = theta, 
      alpha = alpha,
      sd = 2,
      sim_num = b
    )
  }
))

boot_samp <- sampling_results$samples

cl <<- ParallelLogger::makeCluster(ncores)
parallel::clusterExport(cl, varlist = c(
  'ipak',
  'n_repeat',
  'extract_number',
  'create_ar_matrix',
  'boot_samp',
  'ar',
  'theta_of_alpha',
  'split_data',
  'estimate_gee',
  'compute_mu_alpha_jacobian',
  'compute_gee_variance'
))
parallel::clusterEvalQ(cl, ipak('geepack', 'data.table', 'dplyr', 'numDeriv', 'stringr'))


user_results <- ParallelLogger::clusterApply(
  cl,
  1:B,
  function(b){
    data_b <- split_data(boot_samp, b)
    estimate <- estimate_gee(data_b, ar, verbose = F, tol = 1e-06)
    vcov <- compute_gee_variance(estimate, data_b, ar)
    
    return(list(
      theta = estimate$theta,
      alpha = estimate$alpha,
      convergence = estimate$converg,
      vcov = vcov))
  },
  progressBar = TRUE)


gee_results <- ParallelLogger::clusterApply(
  cl,
  1:B,
  function(b){
    geeglm_mod <- geeglm(
        formula = response_value ~ 0 + factor(response_index) + is_sick:factor(response_index),
        data = boot_samp[simulation == b],
        family = gaussian(), id = id, corstr = 'ar'
      )
    
    estimate <- coef(geeglm_mod)
    vcov <- vcov(geeglm_mod)
    alpha_index <- (n_repeat + 1):(2*n_repeat)
    
    return(list(
      theta = estimate[1:n_repeat],
      alpha = estimate[alpha_index],
      convergence = NA,
      vcov = vcov[alpha_index, alpha_index]
      ))
  },
  progressBar = TRUE)

ParallelLogger::stopCluster(cl)

gee_results_agg <- extact_results(gee_results)
user_results_agg <- extact_results(user_results)

##### results #####

estimates <- data.table(
  user_emp = colMeans(user_results_agg$alpha_mat),
  gee_emp = colMeans(gee_results_agg$alpha_mat),
  real = alpha
)


sd_estimates <- data.table(
  user_emp = sqrt_diag(var(user_results_agg$alpha_mat)),
  user = colMeans(user_results_agg$sd_mat),
  gee = colMeans(gee_results_agg$sd_mat),
  gee_emp = sqrt_diag(var(gee_results_agg$alpha_mat))
)

print(estimates)
print(sd_estimates)

lims <- range(sd_estimates)

ggplot(sd_estimates, aes(x = gee, y = user)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  xlim(lims) + ylim(lims)

  

if(do_save){
  file_name <- paste0('main_work/data/enviroments/test_gee_', format(Sys.time(), '%Y-%m-%d-%H%M'), '.RData')
  save.image(file_name)
  message(paste0('file saved: ', file_name))
}
