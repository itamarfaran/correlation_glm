source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$multiplicative_identity
p <- 40
T_thresh <- Tlength <- 115

ARMAdetails <- list(
  ARsick = NULL, MAsick = NULL,
  ARhealth = NULL, MAhealth = NULL
  )
sapply(ARMAdetails, checkInv)

sampleData <- createSamples(
  nH = 30, nS = 30, p = p, Tlength = 115, dim_alpha = 1,
  percent_alpha = 0.3, range_alpha = c(1, 1),
  ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
  ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
  ncores = ncores
  )

# sampleData <- prepare_corrmat_data(
#   link = "main_work/Data/Amnesia_all_AAL.mat",
#   corr_matrix_name = 'corrmats',
#   healthy_index_name = 'CONTROLS',
#   sick_index_name = 'TGA',
#   subset = 1:p
# )
# 
# sampleData <- prepare_corrmat_data(
#   link = "main_work/Data/ADNI_data_AD_CN.mat",
#   corr_matrix_name = 'all.corrmats',
#   healthy_index_name = 'CONTROLS',
#   sick_index_name = 'AD',
#   subset = 1:p
# )
# 
# sampleData <- prepare_corrmat_data(
#   link = "main_work/Data/NMDA_all_data_AAL90.mat",
#   corr_matrix_name = 'group.all',
#   healthy_index_name = 'CONTROLS',
#   sick_index_name = 'NMDA',
#   subset = 1:p
# )

test_corr_mat(sampleData)

Pelet_Cov <- estimateAlpha(
  healthy.data = sampleData$samples$healthy, sick.data = sampleData$samples$sick,
  dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
  T_thresh = T_thresh, updateU = 1, progress = T, linkFun = linkFun)

gee_var <- compute_gee_variance(
  CovObj = Pelet_Cov, sick.data = sampleData$samples$sick,
  linkFun = linkFun, est_mu = FALSE
)
mle_var <- compute_sandwhich_fisher_variance(
  CovObj = Pelet_Cov, sick.data = sampleData$samples$sick,
  linkFun = linkFun, dim_alpha = 1
)

HypTestResMLE <- build_hyp_test(Pelet_Cov, mle_var, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResGEE <- build_hyp_test(Pelet_Cov, gee_var, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
gc()

Pelet_Cov$returns
Pelet_Cov$convergence
c("Est_DF" = Pelet_Cov$Est_N, "Real_DF" = Tlength)
HypTestResMLE$Results
HypTestResGEE$Results

Pelet_Cov_jacknife <- estimateAlpha_jacknife(
  healthy.data = sampleData$samples$healthy, sick.data = sampleData$samples$sick,
  dim_alpha = 1, reg_lambda = 0, var_weights = c(1, 0, 0),
  T_thresh = Tlength, updateU = 1, progress = T, linkFun = linkFun, jack_healthy = TRUE)

alpha_jk_estimate <- colMeans(Pelet_Cov_jacknife$alpha)
alpha_jk_variance <- var(Pelet_Cov_jacknife$alpha)*nrow(Pelet_Cov_jacknife$alpha - 1)
alpha_jk_sds <- sqrt(diag(alpha_jk_variance))
alpha_jk_z <- (alpha_jk_estimate - linkFun$NULL_VAL)/alpha_jk_sds
alpha_jk_pval <- 2*pnorm(abs(alpha_jk_z), lower.tail = FALSE)

pval_plot <-
  data.frame(
    MLE = HypTestResMLE$Results$`P-val`,
    GEE = HypTestResGEE$Results$`P-val`,
    JK = alpha_jk_pval
  ) %>%
  gather(
    key = Method,
    value = pval
  ) %>%
  ggplot(
    aes(
      x = pval,
      fill = Method)
    ) +
  geom_histogram(
    col = 'white',
    binwidth = 0.1,
    boundary = 0,
    show.legend = FALSE 
    ) + 
  geom_hline(yintercept = 0) + 
  geom_hline(
    yintercept = length(alpha_jk_pval)/10,
    col = 'darkgrey',
    size = 1,
    linetype = 2
    ) + 
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    limits = 0:1
    ) + 
  facet_wrap(.~Method) + 
  labs(
    title = 'P-Values with Differnt Variance Estimations',
    x = 'P-Value',
    y = 'Frequency'
    )
pval_plot
save.image('main_work/data/enviroments/full_run_jacknife_simulated_data_multiplicative_link.RData')

