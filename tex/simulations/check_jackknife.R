source('tex/simulations/aux_.R')

infer_jacknife_testing <- function(results, version_){
  if(class(version_) == 'numeric') version_ == as.character(version_)
  sick_ind <- as.logical(results$is_diagnosed)
  
  n_s <- sum(sick_ind)
  n_h <- sum(!sick_ind)
  
  estimate_d <- colMeans(results$alpha[sick_ind,])
  const_d <- (n_s - 1)^2  # /n_s
  var_d <- var(results$alpha[sick_ind,])*const_d
  
  estimate_h <- colMeans(results$alpha[!sick_ind,])
  const_h <- (n_h - 1)^2  # /n_h
  var_h <- var(results$alpha[!sick_ind,])*const_h
  
  estimate <- colMeans(results$alpha)
  # estimate <- (estimate_d + estimate_h)/2
  # var_out <- var_d + var_h
  var_out <- switch(
    version_,
    '1' = (var_d + var_h)/(n_s + n_h),
    '2' = (var_d + var_h)*(1/n_s + 1/n_h),
    '3' = var_d/n_s + var_h/n_h,
    NA) 
    
  
  return(list(
    estimate = estimate,
    variance = var_out
  ))
}


n_sim <- 2*ncores
n <- 50
p <- 24
range_alpha <- c(0.7, 1.1)
percent_alpha <- 0.8
ARMA <- 0.5

n_h <- ceiling(0.5*n)
n_s <- n - n_h
samples <- create_samples(
  n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, 
  percent_alpha = percent_alpha, range_alpha = range_alpha, 
  ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA,
  ncores = ncores
)

estimates <- do.call(rbind, pbmclapply(1:n_sim, function(i){
  out <- estimate_model(
    control_arr = samples$samples[[i]]$healthy,
    diagnosed_arr = samples$samples[[i]]$sick,
    verbose = FALSE)
  return(as.vector(out$alpha))
}, mc.cores = ncores))

jackknife_estimates <- estimate_model_jacknife(
  control_arr = samples$samples[[1]]$healthy,
  diagnosed_arr = samples$samples[[1]]$sick,
  verbose = TRUE, ncores = ncores
)

jackknife_inference_1 <- infer_jacknife_testing(jackknife_estimates, 1)
jackknife_inference_2 <- infer_jacknife_testing(jackknife_estimates, 2)
jackknife_inference_3 <- infer_jacknife_testing(jackknife_estimates, 3)


sds_dt <- data.table(
  jackknife_1 = sqrt_diag(jackknife_inference_1$variance),
  jackknife_2 = sqrt_diag(jackknife_inference_2$variance),
  jackknife_3 = sqrt_diag(jackknife_inference_3$variance),
  empiric = sqrt_diag(var(estimates))
)

sds_dt_long <- melt(sds_dt, id.vars = 'empiric', variable.name = 'method', value.name = 'jackknife')

plt <- ggplot(sds_dt_long, aes(x = empiric, y = jackknife)) + 
  geom_point(shape = 21, col = 'black', fill = 'darkgrey', alpha = 0.7) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(.~method) + 
  theme_user() + 
  labs(
    title = 'Jackknife Estimate of Variance',
    x = 'Empirical Variance', y = 'Jackknife Estimated Variance'
  )

plt <- square_plot(plt)

save.image(file = 'tex/simulations/jacknife_test.RData')

custom_ggsave('jacknife_test.png', plt, width=.8, height=.8)

