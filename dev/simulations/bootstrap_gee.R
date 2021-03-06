source("lib/utils/general_functions.R")
source("lib/utils/simulation_functions.R")
source("lib/model/estimation_functions.R")
source("lib/model/inference_functions.R")

linkFun <- linkFunctions$multiplicative_identity
B <- 5000
ARMAdetails <- list(
  ARsick = NULL, MAsick = NULL,
  ARhealth = NULL, MAhealth = NULL
)

B_uniq <- floor(B/20) # B
combinations2boot <- data.table(
  p = 10*round((40*rbeta(B_uniq, 1, 1.8) + 20)/10),
  n = 10*round(runif(B_uniq, 80, 120)/10),
  sick_obs_percentage = 10*round((0.6*rbeta(B_uniq, 1.5, 1.5) + 0.2)/10, 2),
  Tlength = 10*round((40*rbeta(B_uniq, 1.5, 1.5) + 80)/10),
  ar = round(c(
    runif(B_uniq/4, -0.7, -0.3),
    rep(0, B_uniq/2),
    runif(B_uniq/4, 0.3, 0.7)
  ), 1)
)
combinations2boot <- combinations2boot[rep(1:.N, each = 20)] # commentize
setorder(combinations2boot, p, n, Tlength, ar)

sapply(combinations2boot, function(x) rbind(
  sort(unique(x)),
  tabulate(sort(factor(x)))
))

combinations2boot[,`:=`(
  nh = ceiling((1 - sick_obs_percentage) * n),
  ns = ceiling(sick_obs_percentage * n)
)]

seed <- c(3584, 5345, 3848)

samples <- 
  pbmclapply(
    1:combinations2boot[,.N], function(i){
      samp_ <- create_samples(
        n_sim = 1, #seed = seed,
        n_h = combinations2boot[i, nh],
        n_s = combinations2boot[i, ns],
        p = combinations2boot[i, p],
        Tlength = combinations2boot[i, Tlength],
        dim_alpha = 1, percent_alpha = 0.3, range_alpha = c(1, 1),
        ARsick = combinations2boot[i, ar], MAsick = NULL,
        ARhealth = combinations2boot[i, ar], MAhealth = NULL,
        ncores = 1
      )
      covobj <- estimate_alpha(
        healthy_dt = samp_$samples$healthy,
        sick_dt = samp_$samples$sick,
        dim_alpha = 1, verbose = F, linkFun = linkFun)
      covobj$gee_var_new <- compute_gee_variance(
        covobj, 
        healthy_dt = samp_$samples$healthy,
        sick_dt = samp_$samples$sick
        )
      covobj$gee_var_old <- NA # compute_gee_variance_nosick(covobj, samp_$samples$sick, linkFun = linkFun)
      covobj$mle_var <- NA # compute_sandwhich_fisher_variance(covobj, samp_$samples$sick, linkFun = linkFun)
      covobj$real_theta <- samp_$real_theta
      covobj$real_alpha <- samp_$alpha
      covobj$healthy_data <- convert_corr_array_to_data_matrix_test(samp_$samples$healthy)
      covobj$sick_data <- convert_corr_array_to_data_matrix_test(samp_$samples$sick)
      covobj$Steps <- covobj$Log_Optim <- NULL
      return(covobj)
    }, mc.cores = ncores
  )


combinations2boot[,`:=`(
  actual_sd = sapply(transpose(samples)$alpha, sd_known_mu, mu = 1),
  gee_sd = sapply(transpose(samples)$gee_var_new, mean_sqrt_diag),
  gee_old_sd = NA, # sapply(transpose(samples)$gee_var_old, mean_sqrt_diag),
  mle_sd = NA, # sapply(transpose(samples)$mle_var, mean_sqrt_diag),
  rejected = sapply(1:.N, function(i, sig_level){
    zval <- with(samples[[i]], (alpha - linkFun$NULL_VAL)/sqrt_diag(gee_var_new))
    pval <- 2*pnorm(abs(zval), lower.tail = FALSE)
    return(sum(pval < sig_level))
  }, sig_level = 0.05
  )
)]
combinations2boot[,`:=`(
  type1error = round(rejected/p, 3),
  nh = NULL, ns = NULL
)]

fwrite(combinations2boot, 'main_work/code/gee-bootstrap-app-old/gee_data.csv')
save.image(paste0('main_work/data/enviroments/p_n_bootstrap', format(Sys.time(), '%Y%m%d_%H%M'), '.RData'))


#####

combinations2boot <- fread('main_work/code/gee-bootstrap-app/gee_data.csv')

with(combinations2boot, plot(actual_sd, gee_sd))
with(combinations2boot, cor(actual_sd, gee_sd))
with(combinations2boot, cor(actual_sd, gee_sd, method = 'spearman'))
summary(with(combinations2boot, lm(gee_sd ~ 0 + actual_sd)))
abline(a = 0, b = 1)
hist(combinations2boot$type1error)
summary(combinations2boot$type1error)

cols <- c('gee_sd', 'gee_old_sd', 'mle_sd')
combinations2boot_ratio_long <- 
  combinations2boot[,.(
    p, n, sick_obs_percentage, Tlength,
    actual_sd = 0,
    gee_sd = gee_sd/actual_sd,
    gee_old_sd = gee_old_sd/actual_sd,
    mle_sd = mle_sd/actual_sd
    )] %>% melt(
      measure.vars = (colnames(combinations2boot) %>% (function(x) x[str_detect(x, '_sd')]))
      )
combinations2boot_ratio_long[,.(
  mean = mean(value)#,
  # q10 = quantile(value, 0.1),
  # q90 = quantile(value, 0.9)
  ), by = .(by_clause = p, variable)] %>%
  dcast(by_clause ~ variable) %>%
  (function(dt) dt[order(gee_sd)])

combinations2boot %>%
  mutate(
    m = 0.5*p*(p - 1),
    diff_np = n - p,
    rat_np = n/p,
    diff_nm = n - m,
    rat_nm = n/m,
    gee_sd_ratio = gee_sd/actual_sd
    ) %>% 
  select(
    p, n, 
    sick_obs_percentage,
    Tlength,
    # diff_np, rat_np,
    m, diff_nm, rat_nm,
    actual_sd,
    gee_sd,
    gee_sd_ratio
    ) %>%
  (function(dt, method = 'spearman'){
    m1 <- cor(dt, method = method)
    m2 <- cor.mtest(dt, conf.level = .95, method = method)
    
    return(
      corrplot(
        corr = m1,
        p.mat = m2$p,
        type = "upper",
        method = 'color',
        tl.col = 'black',
        insig = 'blank',
        addCoef.col = 'black',
        number.cex = 0.6,
        addCoefasPercent = FALSE,
        diag = TRUE
      )
    )
  })

