source("Main Work/Code/generalFunctions.R")
source("Main Work/code/estimationFunctions2.R")
source("Main Work/code/simulationFunctions.R")
source("Main Work/code/main_code_noPlots.R")

tt <- rep(Sys.time(), 2)
ncores <- detectCores() - 1 #### Enter here!
if(ncores > 1) requiredFunction <- c("Estimate.Loop", "Estimate.Loop2",
                                     "cor.matrix_to_norm.matrix", "triangle_to_vector","vector_to_triangle",
                                     "create_alpha_mat", "clean_sick", "vnorm", "compute_estimated_N","vector_var_matrix_calc_COR",
                                     "minusloglik", "bootstrapFunction")

Tlength <- 200
ARMAdetails <- list(ARsick = c(0.4, -0.2), MAsick = c(0.4),
                    ARhealth = c(0.2, -0.1), MAhealth = c(0.4))
all(sapply(ARMAdetails, checkInv))

B <- 100
p <- 10
sampleDataB <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                             percent_alpha = 0.4, range_alpha = c(0.6, 0.8),
                             ARsick = ARMAdetails$ARsick, ARhealth = ARMAdetails$ARhealth,
                             MAsick = ARMAdetails$MAsick, MAhealth = ARMAdetails$MAhealth)

bootstrapFunction <- function(b){
  res_unspecified <- Estimate.Loop(Healthy_List = sampleDataB$samples[[b]]$healthy,
                                   Sick_List = sampleDataB$samples[[b]]$sick)
  
  res_specified <- Estimate.Loop2(theta0 = res_unspecified$theta,
                                  alpha0 = res_unspecified$alpha,
                                  healthy.data = sampleDataB$samples[[b]]$healthy,
                                  sick.data = sampleDataB$samples[[b]]$sick,
                                  T_thresh = 10^4, progress = FALSE)
  return(res_specified)
}

tt1 <- Sys.time()
if(ncores > 1){
  buildCL(ncores, c("dplyr", "matrixcalc"), requiredFunction)
  clusterExport(cl = cl, "sampleDataB")
  simuldat <- parLapply(cl = cl, 1:B, bootstrapFunction)
  terminateCL()
} else {
  simuldat <- lapply(1:B, bootstrapFunction)
}
tt1 <- Sys.time() - tt1
alpha_simul <- matrix(nrow = B, ncol = p)
estN_all <- numeric(B)

for(b in 1:B){
  alpha_simul[b,] <- simuldat[[b]]$alpha
  estN_all[b] <- simuldat[[b]]$Est_N
} 
VarAlphaByHess <- ComputeFisher(simuldat[[1]], sampleDataB$samples[[1]]$sick, "Hess") %>% solve
VarAlphaByGrad <- ComputeFisher(simuldat[[1]], sampleDataB$samples[[1]]$sick, "Grad") %>% solve
VarAlphaCombined <- VarAlphaByHess %*% solve(VarAlphaByGrad) %*% VarAlphaByHess

Emp_vs_Theo <- data.frame(TheoreticHess = sqrt(diag(VarAlphaByHess)),
                          TheoreticGrad = sqrt(diag(VarAlphaByGrad)),
                          TheoreticCombined = sqrt(diag(VarAlphaByHess %*% solve(VarAlphaByGrad) %*% VarAlphaByHess)),
                          Empiric = sapply(1:p, function(i) sd(alpha_simul[,i]))) %>%
  mutate(QuotentHess = TheoreticHess/Empiric,
         QuotentGrad = TheoreticGrad/Empiric,
         QuotentCombined = TheoreticCombined/Empiric,
         QuotentBetween = TheoreticHess/TheoreticGrad)

SDErrorByFunction <- Emp_vs_Theo %>% gather(key = Type, value = Thoeretic, -Empiric, -starts_with("Quotent")) %>%
  ggplot(aes(x = Thoeretic, y = Empiric, col = Type)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(intercept = 0, slope = 1, size = 0.9, linetype = 3) + 
  geom_point() + geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, linetype = 2)

BiasDiff <- ggplot(data.frame(Bias = estN_all - Tlength), aes(x = Bias)) +
  geom_histogram(bins = sqrt(B), col = "white", fill = "lightblue") + labs(title = "Bias of Estimated N")
BiasRatio <- ggplot(data.frame(Bias = estN_all/Tlength - 1), aes(x = Bias)) +
  geom_histogram(bins = sqrt(B), col = "white", fill = "lightblue") + labs(title = "Bias of Estimated N")

link3 <- gsub(":", "-", paste0("Main Work/Data/Enviroments/", "fullRunYesARMA ", Sys.time(), ".RData") )

save.image(file = link3)
