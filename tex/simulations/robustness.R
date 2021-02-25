source('tex/simulations/aux_.R')

n <- 50
p <- 24
ARMA <- 0.5
percent_alpha <- 0.8
range_alpha <- c(0.6, 0.95)

create_variance_estimates_different_link <- function(n_sim, simulate_link, name_sim, estimate_link, name_est){
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(0.5*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha, linkFun = simulate_link,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA)
  results <- lapply(
    1:n_sim, function(i) estimate_model(
      control_arr = samples$samples[[i]]$healthy,
      diagnosed_arr = samples$samples[[i]]$sick,
      LinkFunc = estimate_link, verbose = FALSE)
  )
  
  gee_vars <- sapply(1:n_sim, function(i) compute_gee_variance(
    mod = results[[i]],
    control_arr = samples$samples[[i]]$healthy,
    diagnosed_arr = samples$samples[[i]]$sick
  ), simplify = 'array')
  
  sds <- lapply(1:n_sim, function(i) sqrt_diag(gee_vars[,,i]))
  
  
  out <- data.table(
    sims_num = rep(1:n_sim, each = p),
    real = as.vector(samples$alpha),
    estimate = as.vector(sapply(transpose(results)$alpha, as.vector)),
    sd = do.call(c, sds)
  )
  
  out[,`:=`(
    simulate_link = name_sim, 
    estimae_link = name_est,
    z_value = (estimate - estimate_link$null_value)/sd,
    is_null = real == simulate_link$null_value
    )]
  
  return(out)
}

res <- create_variance_estimates_different_link(
  20,
  LinkFunctions$multiplicative_identity, 'additive',
  LinkFunctions$additive_quotent, 'quotent'
  )

res[,`:=`(
  p_value = 2*pnorm(abs(z_value), lower.tail = F),
  is_null_char = ifelse(is_null, 'Null Value', 'Non Null Value')
)]

r_squared <- round(summary(lm(estimate ~ 0 + factor(real), res))$r.squared, 3)

estimates_plot <- ggplot(res, aes(x = real, y = estimate)) + 
  geom_point(shape = 17, position = 'jitter') + 
  annotate('label', x = min(res$real)*1.1, y = 0, label = TeX(paste0('$R^2$: ', r_squared), 'expression')) + 
  labs(
    # title = 'Bias of Alpha Estimate',
    # subtitle = 'Showing Box Plot for Null Values (\u03B1 = 1)',
    x = 'Real Value (Multiplicative)',
    y = 'Estimate (Quotent)'
    ) + 
  theme_bw()

inference_plot <- ggplot(res, aes(x = z_value, y = ..density.., fill = is_null_char)) + 
  geom_histogram(bins = 8*sqrt(length(res)), alpha = .5, position = 'identity', col = 'white') +
  # facet_grid(is_null_char~.) + 
  geom_hline(yintercept = 0) + 
  labs(
    x = 'Z scores \n of misspecified link function',
    y = 'Density',
    fill = ''
  ) +
  theme_user() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(.8,.75),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) + 
  scale_fill_manual(
    values = c('Null Value' = '#A6ACAF', 'Non Null Value' = '#212F3D')
  )

out <- arrangeGrob(estimates_plot, inference_plot, nrow = 2)

custom_ggsave('robustness_link.png', out, width = 2, height = 1.1)
