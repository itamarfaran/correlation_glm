source('tex/simulations/aux_.R')

n_sim <- 2*ncores
p <- 32
DO_CORRECT = TRUE

create_power_comparison <- function(
  sim, n_sim, n, p, percent_alpha, range_alpha, ARMA = 0,
  method = 'BH', sig_level = .05, linkFun = linkFunctions$multiplicative_identity, ncores = 1){
  
  fisher_z <- function(x) 0.5*log((1 + x)/(1 - x))
  
  case <- if(percent_alpha == 0) 'No Effect' else "Effect"
  autocorrelated <- if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(0.5*n)
  n_h <- n - n_s
  
  analyze_sample <- function(i){
    sample_ <- create_samples(n_sim = 1, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                              percent_alpha = percent_alpha, range_alpha = range_alpha,
                              ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA,
                              linkFun = linkFun, enforce_min_alpha = TRUE, ncores = ncores)
    
    results <- estimate_alpha(
        healthy_dt = sample_$samples$healthy,
        sick_dt = sample_$samples$sick,
        linkFun = linkFun, verbose = FALSE)
    
    gee_vars <- compute_gee_variance(
      cov_obj = results,
      healthy_dt = sample_$samples$healthy,
      sick_dt = sample_$samples$sick)
    
    if(DO_CORRECT){
      z_ <- (as.vector(results$alpha) - 1)/sqrt_diag(gee_vars * 1.1)
    } else{
      z_ <- (as.vector(results$alpha) - 1)/sqrt_diag(gee_vars)
    }
    p_ <- 2*pnorm(abs(z_), lower.tail = F)
    p_adj <- p.adjust(p_, method)
    
    to_reject_gee <- sample_$alpha != linkFun$NULL_VAL
    power_gee <- p_adj[to_reject_gee] < sig_level
    error_gee <- p_adj[!to_reject_gee] < sig_level
    
    t_test_mat <- matrix(0, p, p)
    for(i_ in 1:(p-1)) for(j_ in (i_+1):p)
      t_test_mat[i_,j_] <- with(
        sample_$samples,
        t.test(fisher_z(healthy[i_,j_,]), fisher_z(sick[i_,j_,]))$p.value
      )
    t_test_mat[upper.tri(t_test_mat)] <- p.adjust(t_test_mat[upper.tri(t_test_mat)], method)
    t_test_mat <- t_test_mat + t(t_test_mat)
    diag(t_test_mat) <- 1
    
    pvals_t <- t_test_mat[lower.tri(t_test_mat)]
    to_reject_t <- with(sample_, linkFun$FUN(triangle2vector(real_theta), alpha, d=1) != real_theta)
    to_reject_t <- to_reject_t[lower.tri(to_reject_t)]
    
    power_t <- pvals_t[to_reject_t] < sig_level
    error_t <- pvals_t[!to_reject_t] < sig_level
    
    out <- data.table(
      sim_num = i,
      t_true_null = sum(!to_reject_t),
      t_fp = sum(error_t),
      t_true_alt = sum(to_reject_t),
      t_tp = sum(power_t),
      gee_true_null = sum(!to_reject_gee),
      gee_fp = sum(error_gee),
      gee_true_alt = sum(to_reject_gee),
      gee_tp = sum(power_gee)
    )
    
    return(out)
  }
  
  out <- do.call(rbind, mclapply(1:n_sim, analyze_sample, mc.cores = ncores))
  out[,`:=`(
    sim = sim, n = n, p = p, percent_alpha = percent_alpha, min_alpha = min(range_alpha),
    autocorrelated = autocorrelated, case = case
  )]
  
  return(out)
}


min_alpha <- c(0.8, 0.85, 0.9, 0.95) 
examples <- unique(
  rbind(
    expand.grid(
      min_alpha = min_alpha,
      n = 100,
      percent_alpha = c(.05, .2)
    ),
    expand.grid(
      min_alpha = min_alpha,
      n = c(60, 120),
      percent_alpha = .1
    )
  )
)

file_loc <- if(DO_CORRECT) 'tex/simulations/power_t_correct.RData' else 'tex/simulations/power_t.RData'
if(file.exists(file_loc)){
  load(file_loc)
} else {
  out <- pblapply(seq_len(nrow(examples)), function(i) create_power_comparison(
    sim = i, n_sim = n_sim, n = examples[i, 2], p = p, percent_alpha = examples[i, 3],
    range_alpha = c(examples[i, 1], 1), ncores = ncores
    ))
  
  out <- do.call(rbind, out)
  
  save(out, file = file_loc)
}

# out <- out[(percent_alpha %in% c(.05, .2)) | (n %in% c(60, 120))]

out[,`:=`(
  t_fdr = t_fp/pmax(t_fp + t_tp, 1),
  t_power = t_tp/t_true_alt,
  t_1rejected = 1*(t_tp > 0),
  gee_fdr = gee_fp/pmax(gee_fp + gee_tp, 1),
  gee_power = gee_tp/gee_true_alt,
  gee_1rejected = 1*(gee_tp > 0)
)]


id_cols <- c('sim_num', 'n', 'p', 'percent_alpha', 'min_alpha', 'autocorrelated', 'case')
value_cols <- c('t_fdr', 't_power', 't_1rejected', 'gee_fdr', 'gee_power', 'gee_1rejected')
cols <- c(id_cols, value_cols)
toplot <- melt(out[,..cols], id.vars = id_cols)

cols <- c('method', 'rate')
toplot[,(cols) := asplit(do.call(rbind, str_split(toplot[,variable], '_')), 2)]
toplot[,variable := NULL]
toplot[,method := ifelse(method == 't', 'T Test', 'GEE')]

toplot_groupby <- toplot[,.(value = mean(value), lower = quantile(value, .025), upper = quantile(value, .975), .N),
                         by = .(n = factor(n), p, min_alpha, percent_alpha = factor(percent_alpha), method, rate)]
toplot_groupby[,se := sqrt(value*(1 - value)/N)]
toplot_groupby[,`:=`(lower = pmax(value - qnorm(.975) * se, 0), upper = pmin(value + qnorm(.975) * se, 1))]
toplot_groupby[,alpha_inv := 1 - min_alpha]

# todo: minimal decay (1-alpha) in axis
plt1 <- toplot_groupby %>%
  filter(rate %in% c('power', '1rejected'), percent_alpha == .1) %>% 
  mutate(rate = ifelse(rate == 'power', 'Statistical Power', 'Prob. to Reject Global Null')) %>%
  # mutate(value = ifelse(value == 0, 10^-5, value),
  #        lower = ifelse(lower == 0, 10^-5, lower),
  #        upper = ifelse(upper == 0, 10^-5, upper)) %>% 
  ggplot(aes(x = alpha_inv, y = value, ymin = lower, ymax = upper, shape = method, linetype = method, color = n)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c('darkgrey', 'black')) + 
  # geom_errorbar(width = .01) + 
  facet_grid(rate ~ .) + 
  labs(x = TeX('Maximal Decay Effect ($\\max_j \\[ 1-\\alpha_{j}\\]$)'),
       y = 'Rate', col = 'Sample Size', linetype = 'Method', shape = 'Method') +
  theme_user() +
  theme(legend.position = 'bottom', legend.box = 'vertical', legend.margin = margin())


plt2 <- toplot_groupby %>%
  filter(rate %in% c('power', '1rejected'), n == 100) %>% 
  mutate(rate = ifelse(rate == 'power', 'Statistical Power', 'Prob. to Reject Global Null')) %>%
  # mutate(value = ifelse(value == 0, 10^-5, value),
  #        lower = ifelse(lower == 0, 10^-5, lower),
  #        upper = ifelse(upper == 0, 10^-5, upper)) %>% 
  ggplot(aes(x = alpha_inv, y = value, ymin = lower, ymax = upper, shape = method, linetype = method, color = percent_alpha)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c('darkgrey', 'black')) + 
  # geom_errorbar(width = .01) + 
  facet_grid(rate ~ .) + 
  labs(x = TeX('Maximal Decay Effect ($\\max_j \\[ 1-\\alpha_{j}\\]$)'), 
       y = 'Rate', col = TeX('% of Non-Null $\\alpha$-s'), linetype = 'Method', shape = 'Method') +
  theme_user() +
  theme(legend.position = 'bottom', legend.box = 'vertical', legend.margin = margin())


if(DO_CORRECT){
  custom_ggsave('power_t1_correct.png', plt1, width = 1.2, height = 1.2)
  custom_ggsave('power_t2_correct.png', plt2, width = 1.2, height = 1.2)
} else {
  custom_ggsave('power_t1.png', plt1, width = 1.2, height = 1.2)
  custom_ggsave('power_t2.png', plt2, width = 1.2, height = 1.2)
}
