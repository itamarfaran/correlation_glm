source("main_work/code/01_general_functions.R")
source("main_work/code/02_simulation_functions.R")
source("main_work/code/03_estimation_functions.R")
source("main_work/code/04_inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
B <- 500

bootstrap_gee <- function(p, n, percent_alpha = 0.3, range_alpha = c(1, 1),
                          sick_obs_percentage = 0.5, t_length = 115,
                          ar = 0, sig_level = 0.05,
                          sim_num = NULL, ncores = 1, verbose = FALSE){
  lapply_ <- if(verbose) pbmcapply::pbmclapply else parallel::mclapply
  
  sample_ <- create_samples(
    n_sim = 1,
    n_h = round(n*(1 - sick_obs_percentage)),
    n_s = round(n*sick_obs_percentage),
    p = p,
    Tlength = t_length,
    dim_alpha = 1,
    percent_alpha = percent_alpha, range_alpha = range_alpha,
    ARsick = ar, MAsick = NULL,
    ARhealth = ar, MAhealth = NULL,
    linkFun = linkFun,
    ncores = ncores
  )
  
  covobj <- estimate_alpha(
    healthy_dt = sample_$samples$healthy,
    sick_dt = sample_$samples$sick,
    dim_alpha = 1, verbose = F, linkFun = linkFun)
  
  gee_var <- compute_gee_variance(
    covobj, 
    healthy_dt = sample_$samples$healthy,
    sick_dt = sample_$samples$sick
  )
  
  covobj$steps <- covobj$Log_Optim <- NULL
  
  out <- data.table(
    sim_num = sim_num,
    p = p,
    n = n,
    percent_alpha = percent_alpha,
    min_alpha = min(range_alpha),
    sick_obs_percentage = sick_obs_percentage,
    Tlength = t_length,
    ar = ar,
    confidence = 1 - sig_level,
    alpha = as.vector(sample_$alpha),
    alpha_est = as.vector(covobj$alpha),
    gee_sd = sqrt_diag(gee_var)
  )
  return(out)
}

combinations2boot <- data.table(
  p = 10*round((40*rbeta(B, 1, 1.8) + 20)/10),
  n = 10*round(runif(B, 80, 120)/10),
  sick_obs_percentage = 10*round((0.6*rbeta(B, 1.5, 1.5) + 0.2)/10, 2),
  t_length = 10*round((40*rbeta(B, 1.5, 1.5) + 80)/10),
  percent_alpha = sample(c(0, 0.1, 0.2, 0.3), B, TRUE),
  min_alpha = sample(c(0.7, 0.8, 0.9), B, TRUE),
  ar = round(c(
    runif(B/4, -0.7, -0.3),
    rep(0, B/2),
    runif(B/4, 0.3, 0.7)
  ), 1)
)
combinations2boot[percent_alpha == 0, min_alpha := 1]

results <- do.call(rbind, pbmclapply(1:B, function(b){
  with(combinations2boot[b], bootstrap_gee(
    p = p,
    n = n,
    sick_obs_percentage = sick_obs_percentage,
    t_length = t_length,
    percent_alpha = percent_alpha,
    range_alpha = c(min_alpha, 1),
    ar = ar,
    ncores = ncores
  ))
}, mc.cores = 1))

results[, zval := (alpha_est - linkFun$NULL_VAL)/gee_sd]
results[, pval := 2*pnorm(abs(zval), lower.tail = F)]

beta_mle <- function(x, ini = c(1, 1)){
  minloglik <- function(par, x) -sum(dbeta(x, par[1], par[2], log = T))
  out <- optim(ini, minloglik, x = x, method = 'BFGS', hessian = T)
  return(out)
}

beta_hyp_test <- function(x, pars = c(1, 1), sig_level = 0.05){
  mle <- beta_mle(x, pars)

  dist <- mle$par - pars
  chisq <- as.vector(t(dist) %*% mle$hessian %*% dist)
  pval <- pchisq(chisq, 1)
  ci <- mle$par %o% rep(1, 3) + qnorm(1 - sig_level/4)* sqrt(diag(solve(mle$hessian))) %o% (-1:1)
  colnames(ci) <- c(
    paste0('lower', round(sig_level*100), '%'),
    'estimate',
    paste0('upper', round(sig_level*100), '%')
  )
  rownames(ci) <- c('alpha', 'beta')
  return(list(pvalue = pval, ci = ci))
}

# density_dt <- data.table(
#   hyp = rep(c('Null', 'Alt.'), each = 1000),
#   x = rep(seq(0, 1, length.out = 1000), 2)
# )
# 
# null_beta <- beta_mle(results[alpha == 1, pval])
# alt_beta <- beta_mle(results[alpha != 1, pval])
# density_dt[hyp == 'Null', y := dbeta(x, null_beta$par[1], null_beta$par[2])]
# density_dt[hyp == 'Alt.', y := dbeta(x, alt_beta$par[1], alt_beta$par[2])]

plt <- ggplot(results, aes(x = pval, y = ..density..)) + 
  geom_histogram(
    col = 'white',
    fill = 'lightblue',
    bins = 11, #2*log(nrow(results)),
    boundary = 0
    ) +
  geom_hline(yintercept = 0) +
  geom_density(
    size = 1,
    linetype = 2,
    kernel = 'gaussian'
    ) + 
  xlim(0:1) + 
  facet_grid(
    ifelse(alpha == 1, 'Null', 'Alt.')~.,
    scales = 'free_y'
    ) + 
  labs(
    title = 'Distribution of P-values by Hypothesis',
    x = 'P-Values', y = 'Density'
    ) + 
  theme_minimal() + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
    )

ggsave('main_work/code/gee-bootstrap-app/pval_hist.png', plt)

fwrite(results, 'main_work/code/gee-bootstrap-app/gee_pval_boot.csv')
save.image(paste0('main_work/data/enviroments/gee_pval_bootstrap', format(Sys.time(), '%Y%m%d_%H%M'), '.RData'))
