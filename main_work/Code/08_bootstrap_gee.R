source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- c(5, 10) # c(5, 10, 20, 30, 40, 50, 60)
n <- c(50, 80) # c(50, 80, 100, 120)
sick_obs_percentage <- c(0.5, 0.6) # c(0.2, 0.35, 0.5, 0.65, 0.8)
reps <- 3
T_thresh <- Tlength <- 115
ncores <- ifelse(
  .Platform$OS.type == 'windows', 1,
  max(1, detectCores() - 2)
)
seed <- c(3584, 5345, 3848)

ARMAdetails <- list(
  ARsick = NULL, MAsick = NULL,
  ARhealth = NULL, MAhealth = NULL
)

combinations2boot <- data.table(
  expand.grid(p, n, sick_obs_percentage)
)
colnames(combinations2boot) <- c('p', 'n', 'sick_obs_percentage')

# combinations2boot <- data.table(do.call(
#   rbind, replicate(
#     reps,
#     expand.grid(p, n),
#     simplify = FALSE)
#   ))
# colnames(combinations2boot) <- c('p', 'n')

combinations2boot[,`:=`(
  nh = ceiling((1 - sick_obs_percentage) * n),
  ns = ceiling(sick_obs_percentage * n),
  index = 1:.N,
  rep = ceiling(100/p)
)]
setorder(combinations2boot, p, n)
combinations2boot <- combinations2boot[rep(index, rep)]
combinations2boot[,`:=`(
  index = NULL,
  rep = NULL
)]

samples <- 
  pbmclapply(
    1:combinations2boot[,.N],
    function(i){
      samp_ <- createSamples(
        B = 1, #seed = seed,
        nH = combinations2boot[i, nh],
        nS = combinations2boot[i, ns],
        p = combinations2boot[i, p],
        Tlength = Tlength, dim_alpha = 1,
        percent_alpha = 0.3, range_alpha = c(1, 1),
        ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
        ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
        ncores = 1
      )
      covobj <- estimateAlpha(
        healthy.data = samp_$samples$healthy,
        sick.data = samp_$samples$sick,
        dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
        T_thresh = T_thresh, updateU = 1, progress = F, linkFun = linkFun)
      covobj$gee_var <- compute_gee_variance(covobj, samp_$samples)
      covobj$real.theta <- samp_$real.theta
      covobj$real.alpha <- samp_$alpha
      covobj$healthy_data <- corr_mat_array2normal_data_mat(samp_$samples$healthy)
      covobj$sick_data <- corr_mat_array2normal_data_mat(samp_$samples$sick)
      covobj$Steps <- covobj$Log_Optim <- NULL
      return(covobj)
      },
    mc.cores = ncores
  )


combinations2boot[,`:=`(
  nh = NULL,
  ns = NULL,
  actual_sd = sapply(transpose(samples)$alpha, sd),
  mean_est_sd = sapply(transpose(samples)$gee_var, function(x) mean(sqrt(diag(x)))),
  rejected = sapply(1:.N, function(i, sig_level){
    zval <- with(samples[[i]], (alpha - linkFun$NULL_VAL)/sqrt(diag(gee_var)))
    pval <- 2*pnorm(abs(zval), lower.tail = FALSE)
    return(sum(pval < sig_level))
    }, sig_level = 0.05
  )
)]
combinations2boot[,actual_est_ratio := actual_sd/mean_est_sd]

ggplot(combinations2boot, aes(
  x = p, y = actual_est_ratio,
  col = factor(n))#, shape = factor(sick_obs_percentage))
  ) +
  geom_point() + 
  geom_hline(yintercept = 0, size = 1) + 
  geom_hline(yintercept = 1, size = 1, col = 'darkgrey', linetype = 2) + 
  stat_summary(fun.y = mean, geom = "line", size = 1.2, linetype = 1) + 
  labs(
    x = 'P',
    y = 'Ratio between Estimated and Actual SD',
    col = 'N',
    shape = 'Sick Obs / Total Obs'
  )

lubridate::round_date(Sys.time(), 'minute')
save.image(paste0('main_work/data/enviroments/p_n_bootstrap', format(Sys.time(), '%Y%m%d_%H%M'), '.RData'))