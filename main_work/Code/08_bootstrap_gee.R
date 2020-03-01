source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/041_inferenceFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$multiplicative_identity
B <- 2000
ARMAdetails <- list(
  ARsick = NULL, MAsick = NULL,
  ARhealth = NULL, MAhealth = NULL
)

combinations2boot <- data.table(
  # p = 5*round((50*rbeta(B, 1, 3) + 10)/5),
  # p = rep(c(20, 40), each = B/2),
  p = 5*round(runif(B, 20, 40)/5),
  n = 10*round(runif(B, 50, 200)/10),
  sick_obs_percentage = 5*round((0.7*rbeta(B, 3, 3) + 0.1)/5, 2),
  Tlength = 5*round((70*rbeta(B, 3, 3) + 80)/5)
)
setorder(combinations2boot, p, n)

combinations2boot[,`:=`(
  nh = ceiling((1 - sick_obs_percentage) * n),
  ns = ceiling(sick_obs_percentage * n)
)]

seed <- c(3584, 5345, 3848)

samples <- 
  pbmclapply(
    1:combinations2boot[,.N],
    function(i){
      samp_ <- createSamples(
        B = 1, #seed = seed,
        nH = combinations2boot[i, nh],
        nS = combinations2boot[i, ns],
        p = combinations2boot[i, p],
        Tlength = combinations2boot[i, Tlength],
        dim_alpha = 1, percent_alpha = 0.3, range_alpha = c(1, 1),
        ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
        ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
        ncores = 1
      )
      covobj <- estimateAlpha(
        healthy.data = samp_$samples$healthy,
        sick.data = samp_$samples$sick,
        dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
        T_thresh = combinations2boot[i, Tlength], updateU = 1, progress = F, linkFun = linkFun)
      covobj$gee_var_new <- compute_gee_variance(covobj, samp_$samples)
      covobj$gee_var_old <- compute_gee_variance_nosick(covobj, samp_$samples$sick, linkFun = linkFun)
      covobj$mle_var <- compute_sandwhich_fisher_variance(covobj, samp_$samples$sick, linkFun = linkFun)
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
  actual_sd = sapply(transpose(samples)$alpha, sd),
  gee_sd = sapply(transpose(samples)$gee_var_new, mean_sqrt_diag),
  gee_old_sd = sapply(transpose(samples)$gee_var_old, mean_sqrt_diag),
  mle_sd = sapply(transpose(samples)$mle_var, mean_sqrt_diag),
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

fwrite(combinations2boot, 'main_work/Code/gee-bootstrap-app/gee_data.csv')
save.image(paste0('main_work/Data/Enviroments/p_n_bootstrap', format(Sys.time(), '%Y%m%d_%H%M'), '.RData'))
