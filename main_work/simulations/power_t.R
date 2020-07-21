source('main_work/simulations/auxilary_functions.R')

n_sim = ncores
sim = 20
p = 32

create_power_comparison <- function(
  sim, n_sim, n, p, percent_alpha, range_alpha, ARMA = 0,
  method = 'BH', sig_level = .05, fraction = 1/5, ncores = 1){
  fisher_z <- function(x) 0.5*log((1 + x)/(1 - x))
  
  case = if(percent_alpha == 0) 'No Effect' else "Effect"
  autocorrelated = if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(0.5*n)
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
  
  gee_vars <- pbmclapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick
  ), mc.cores = ncores)
  
  compare <- function(i){
    z_ <- (as.vector(results[[i]]$alpha) - 1)/sqrt_diag(gee_vars[[i]])
    p_ <- 2*pnorm(abs(z_), lower.tail = F)
    p_adj <- p.adjust(p_, method)
    out_gee <- p_adj < sig_level
    
    t_test_mat <- matrix(0, p, p)
    for(i_ in 1:(p-1)) for(j_ in (i_+1):p)
      t_test_mat[i_,j_] <- with(
        samples$samples[[i]],
        t.test(fisher_z(healthy[i_,j_,]), fisher_z(sick[i_,j_,]))$p.value
      )
    
    t_test_mat[lower.tri(t_test_mat)] <- p.adjust(t_test_mat[lower.tri(t_test_mat)], method)
    t_test_mat <- t_test_mat + t(t_test_mat)
    diag(t_test_mat) <- 1
    out_t <- colSums(t_test_mat < sig_level) > fraction*p
    
    out <- data.table(sim_num = i, voxel = 1:p, rejected_gee = out_gee, rejected_t = out_t)
    return(out)
  }
  
  out <- do.call(rbind, lapply(1:n_sim, compare))
  out[,`:=`(
    sim = sim, estimated_alpha = do.call(c, transpose(results)$alpha),
    real_alpha = rep(samples$alpha, times = n_sim),
    n = n, p = p, percent_alpha = percent_alpha, min_alpha = min(range_alpha),
    sd = do.call(c, lapply(gee_vars, sqrt_diag)),
    autocorrelated = autocorrelated, case = case
  )]
  
  return(out)
}


examples <- expand.grid(n = c(60, 90, 120), percent_alpha = c(0.1, 0.2, 0.3), min_alpha = c(0.85, 0.9, 0.95))

out <- lapply(1:nrow(examples), function(i)
  lapply(
    1:sim, create_power_comparison,
    n_sim = n_sim, n = examples[i, 1], p = p, percent_alpha = examples[i, 2], range_alpha = c(examples[i, 3], 1),
    ncores = ncores
    )
  )

out <- do.call(rbind, lapply(out, do.call, what=rbind))

save(out, file = 'main_work/simulations/power_t.RData')


id.vars <- names(out)[!startsWith(names(out), 'rejected')]
out2 <- melt(out, id.vars=id.vars, variable.name = 'Method')
out2[,Method := ifelse(Method == 'rejected_gee', 'GEE', 'T Test')]

add_same_stuff <- function(plt){
  plt <- plt + 
    scale_y_log10() + 
    # geom_hline(yintercept = 0) + 
    scale_fill_manual(values = c('GEE' = '#505050', 'T Test' = '#DCDCDC')) + 
    theme_user() + theme(
      legend.position = 'none',
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  return(plt)
}

p1 <- out2[real_alpha != 1, mean(value), by = .(sim, percent_alpha, Method)] %>% 
  ggplot(aes(x = percent_alpha, y = V1, fill = Method, group = interaction(percent_alpha, Method))) + 
  geom_boxplot(position = position_dodge(width = .03), alpha = .6) + 
  labs(title = '% Non-Null Parameters') + scale_x_continuous(labels = scales::percent)
p1 <- add_same_stuff(p1); p1
  
p2 <- out2[real_alpha != 1, mean(value), by = .(sim, min_alpha, Method)] %>% 
  ggplot(aes(x = min_alpha, y = V1, fill = Method, group = interaction(min_alpha, Method))) + 
  geom_boxplot(position = position_dodge(width = .02), alpha = .6) + 
  labs(title = 'Minimal Value of Alpha')
p2 <- add_same_stuff(p2); p2

p3 <- out2[real_alpha != 1, mean(value), by = .(sim, n, Method)] %>% 
  ggplot(aes(x = n, y = V1, fill = Method, group = interaction(n, Method))) + 
  geom_boxplot(position = position_dodge(width = 9), alpha = .6) + 
  labs(title = '# Subjects')
p3 <- add_same_stuff(p3); p3

out3 <- arrangeGrob(p1, p2, p3, nrow=3)
plot(out3)

custom_ggsave('power_r.png', out3)