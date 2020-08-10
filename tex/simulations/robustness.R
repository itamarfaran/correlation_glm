source('tex/simulations/aux.R')

n = 50
p = 12
ARMA = 0.5
percent_alpha = 0.4
range_alpha = c(0.6, 1)

create_variance_estimates_different_link <- function(n_sim, simulate_link, name_sim, estimate_link, name_est){
  if (ARMA == 0) ARMA <- NULL
  n_s <- ceiling(0.5*n)
  n_h <- n - n_s
  
  samples <- create_samples(n_sim = n_sim, n_h = n_h, n_s = n_s, p = p, Tlength = 115,
                            percent_alpha = percent_alpha, range_alpha = range_alpha, linkFun = simulate_link,
                            ARsick = ARMA, ARhealth = ARMA, MAsick = ARMA, MAhealth = ARMA)
  results <- lapply(
    1:n_sim, function(i) estimate_alpha(
      healthy_dt = samples$samples[[i]]$healthy,
      sick_dt = samples$samples[[i]]$sick,
      linkFun = estimate_link, verbose = FALSE)
  )
  
  gee_vars <- sapply(1:n_sim, function(i) compute_gee_variance(
    cov_obj = results[[i]],
    healthy_dt = samples$samples[[i]]$healthy,
    sick_dt = samples$samples[[i]]$sick
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
    z_value = (estimate - estimate_link$NULL_VAL)/sd,
    is_null = real == simulate_link$NULL_VAL
    )]
  
  return(out)
}

res <- create_variance_estimates_different_link(
  20,
  linkFunctions$multiplicative_identity, 'additive',
  linkFunctions$additive_quotent, 'quotent'
  )

res[,`:=`(
  p_value = 2*pnorm(abs(z_value), lower.tail = F),
  is_null_char = ifelse(is_null, 'Null Value', 'Non Null Value')
)]

ggplot(res, aes(x = z_value, y = ..density..)) + 
  geom_histogram(bins = 8*sqrt(length(res)), fill = 'lightgrey', col = 'white') +
  facet_grid(is_null_char~.) + 
  geom_hline(yintercept = 0) + 
  theme_user()

ggplot(mapping = aes(x = real, y = estimate)) + 
  geom_point(data = res[is_null == FALSE], shape = 17, position = 'jitter') + 
  geom_boxplot(data = res[is_null == TRUE], width = 0.1, outlier.shape = 1) + 
  ggtitle('Bias of Alpha Estimate', 'Showing Box Plot for Null Values (\u03B1 = 1)') + 
  theme_bw()

