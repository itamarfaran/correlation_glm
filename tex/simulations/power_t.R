source('tex/simulations/aux_.R')

#  todo: plot P(At least one rejected) instead of E[How Much Rejected]

n_sim <- ncores
p <- 32

create_power_comparison <- function(
  sim, n_sim, n, p, percent_alpha, range_alpha, ARMA = 0,
  method = 'BH', sig_level = .05, linkFun = linkFunctions$multiplicative_identity, ncores = 1){
  fisher_z <- function(x) 0.5*log((1 + x)/(1 - x))
  
  case <- if(percent_alpha == 0) 'No Effect' else "Effect"
  autocorrelated <- if(ARMA == 0) 'Not Autocorrelated' else 'Autocorrelated'
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(0.5*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA,
                            linkFun = linkFun, ncores = ncores)
  results <- mclapply(
    1:n_sim, function(i) estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      linkFun = linkFun, verbose = FALSE), mc.cores = ncores
  )
  
  gee_vars <- mclapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick), mc.cores = ncores)
  
  compare <- function(i){
    z_ <- (as.vector(results[[i]]$alpha) - 1)/sqrt_diag(gee_vars[[i]])
    p_ <- 2*pnorm(abs(z_), lower.tail = F)
    p_adj <- p.adjust(p_, method)
    
    to_reject_gee <- samples$alpha != linkFun$NULL_VAL
    power_gee <- p_adj[to_reject_gee] < sig_level
    error_gee <- p_adj[!to_reject_gee] < sig_level
    
    t_test_mat <- matrix(0, p, p)
    for(i_ in 1:(p-1)) for(j_ in (i_+1):p)
      t_test_mat[i_,j_] <- with(
        samples$samples[[i]],
        t.test(fisher_z(healthy[i_,j_,]), fisher_z(sick[i_,j_,]))$p.value
      )
    t_test_mat[upper.tri(t_test_mat)] <- p.adjust(t_test_mat[upper.tri(t_test_mat)], method)
    t_test_mat <- t_test_mat + t(t_test_mat)
    diag(t_test_mat) <- 1
    
    pvals_t <- t_test_mat[lower.tri(t_test_mat)]
    to_reject_t <- with(samples, linkFun$FUN(triangle2vector(real_theta), alpha, d=1) != real_theta)
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
  
  out <- do.call(rbind, lapply(1:n_sim, compare))
  out[,`:=`(
    sim = sim, n = n, p = p, percent_alpha = percent_alpha, min_alpha = min(range_alpha),
    autocorrelated = autocorrelated, case = case
  )]
  
  return(out)
}


examples <- expand.grid(
  n = c(60, 80, 100, 120),
  percent_alpha = c(0.05, 0.1, 0.15, 0.2),
  min_alpha = c(0.8, 0.85, 0.9, 0.95)
  )

file_loc <- 'tex/simulations/power_t.RData'
if(file.exists(file_loc)){
  load(file_loc)
} else {
  out <- pblapply(seq_len(nrow(examples)), function(i) create_power_comparison(
    sim = i, n_sim = n_sim, n = examples[i, 1], p = p, percent_alpha = examples[i, 2],
    range_alpha = c(examples[i, 3], 1), ncores = ncores
    ))
  
  out <- do.call(rbind, out)
  
  save(out, file = file_loc)
}

out[,`:=`(
  t_fdr = t_fp/pmax(t_fp + t_tp, 1),
  t_power = t_tp/t_true_alt,
  gee_fdr = gee_fp/pmax(gee_fp + gee_tp, 1),
  gee_power = gee_tp/gee_true_alt
)]


id_cols <- c('sim_num', 'n', 'p', 'percent_alpha', 'min_alpha', 'autocorrelated', 'case')
value_cols <- c('t_fdr', 't_power', 'gee_fdr', 'gee_power')
cols <- c(id_cols, value_cols)
toplot <- melt(out[,..cols], id.vars = id_cols)

cols <- c('method', 'rate')
toplot[,(cols) := asplit(do.call(rbind, str_split(toplot[,variable], '_')), 2)]
toplot[,variable := NULL]
toplot[,method := ifelse(method == 't', 'T Test', 'GEE')]


plot_by <- function(x, title, scales, width = NULL, type = 'power'){
  toplot %>%
    filter(rate == type) %>% 
    group_by(.dots = c('method', x)) %>% 
    summarise(value = median(value)) %>%
    as.data.table() ->
    labels
  
  interactions_all <- list(
    x = toplot[rate == type, get(x)],
    y = toplot[rate == type, method]
  )
  interactions_lab <- list(
    x = labels[, get(x)],
    y = labels[, method]
  )
  
  p_out <-
    ggplot(toplot[rate == type], aes_string(x = x, y = 'value', fill = 'method')) +
    geom_boxplot(
      aes(group = interaction(interactions_all$x, interactions_all$y)),
      position = position_dodge(width = width), alpha = .6) +
    # geom_label(aes(label = method, group = NULL), labels, position = position_dodge(width = width), fill = 'white')  +
    labs(title = title) +
    scale_x_continuous(labels = scales) + 
    scale_y_log10(limits = c(.001, 1)) +
    # geom_hline(yintercept = 0) +
    scale_fill_manual(values = c('GEE' = '#505050', 'T Test' = '#DCDCDC')) + 
    theme_user() + theme(
      legend.position = 'none',
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )
  return(p_out)
}

out3 <- arrangeGrob(
  plot_by('percent_alpha', '% Non-Null Parameters', scales::percent, width = .04),
  plot_by('min_alpha', 'Minimal Value of Alpha', scales::number, width = .02),
  plot_by('n', '# Subjects', scales::number, width = 9),
  nrow=3)
plot(out3)

custom_ggsave('power_t.png', out3, width = 2, height = 1.5)
