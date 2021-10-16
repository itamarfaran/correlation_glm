source('tex/simulations/aux_.R')

n_repeat <- 5
n_list <- c(40, 80, 120)
p_list <- c(25, 40, 55)

n_sim <- ncores
range_alpha <- c(0.7, 1.1)
percent_alpha <- 0.8
ARMA <- 0.5

results <- list()
i <- 1

for(n in n_list){
  for(p in p_list){
    for(j in n_repeat){
      n_h <- ceiling(0.5*n)
      n_s <- n - n_h
      
      samples <- create_samples(
        n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, 
        percent_alpha = percent_alpha, range_alpha = range_alpha, 
        ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA,
        ncores = ncores
      )
      
      estimates <- do.call(rbind, pbmclapply(1:n_sim, function(k){
        out <- estimate_model(
          control_arr = samples$samples[[k]]$healthy,
          diagnosed_arr = samples$samples[[k]]$sick,
          verbose = FALSE)
        return(as.vector(out$alpha))
      }, mc.cores = ncores))
      
      jackknife_estimates <- estimate_model_jacknife(
        control_arr = samples$samples[[1]]$healthy,
        diagnosed_arr = samples$samples[[1]]$sick,
        verbose = TRUE, ncores = ncores
      )
      
      jackknife_inference <- infer_jackknife(jackknife_estimates)
      
      
      sds_dt <- data.table(
        sim_num = j,
        n = n,
        p = p,
        jackknife = sqrt_diag(jackknife_inference$variance),
        empiric = sqrt_diag(var(estimates))
      )
      
      results[[i]] <- sds_dt
      i <- i + 1
    }
  }
}

sds_dt <- do.call(rbind, results)

plt <- ggplot(sds_dt, aes(x = empiric, y = jackknife)) + 
  geom_point(shape = 21, col = 'black', fill = 'darkgrey', alpha = 0.7) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(p ~ n, labeller = label_both) + 
  theme_user() + 
  labs(
    title = 'Jackknife Estimates for Standard Deviations',
    x = 'Empirical Standard Deviations', y = 'Jackknife Estimated Standard Deviations'
  )

# plt <- square_plot(plt)

save.image(file = 'tex/simulations/jacknife.RData')

custom_ggsave('jacknife.png', plt, width=2, height=1.2)

