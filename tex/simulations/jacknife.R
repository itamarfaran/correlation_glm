source('tex/simulations/aux_.R')

n_repeat <- 5
n_sim <- ncores
n <- 50
p <- 32
range_alpha <- c(0.7, 1.1)
percent_alpha <- 0.8
ARMA <- 0.5

n_h <- ceiling(0.5*n)
n_s <- n - n_h

results <- list()

for(i in seq_len(n_repeat)){
  samples <- create_samples(
    n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, 
    percent_alpha = percent_alpha, range_alpha = range_alpha, 
    ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA,
    ncores = ncores
  )
  
  estimates <- do.call(rbind, pbmclapply(1:n_sim, function(i){
    out <- estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      verbose = FALSE)
    return(as.vector(out$alpha))
  }, mc.cores = ncores))
  
  jackknife_estimates <- estimate_alpha_jacknife(
    healthy_dt = samples$samples[[1]]$healthy,
    sick_dt = samples$samples[[1]]$sick,
    verbose = TRUE, ncores = ncores
  )
  
  jackknife_inference <- infer_jacknife(jackknife_estimates)
  
  
  sds_dt <- data.table(
    sim_num = i,
    jackknife = sqrt_diag(jackknife_inference$variance),
    empiric = sqrt_diag(var(estimates))
  )
  
  results[[i]] <- sds_dt
}

sds_dt <- do.call(rbind, results)

plt <- ggplot(sds_dt, aes(x = empiric, y = jackknife)) + 
  geom_point(shape = 21, col = 'black', fill = 'darkgrey', alpha = 0.7) + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_user() + 
  labs(
    title = 'Jackknife Estimate of Variance',
    x = 'Empirical Variance', y = 'Jackknife Estimated Variance'
  )

plt <- square_plot(plt)

save.image(file = 'tex/simulations/jacknife.RData')

custom_ggsave('jacknife.png', plt, width=.8, height=.8)

