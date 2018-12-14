source("Main Work/Code/generalFunctions.R")
source("Main Work/Code/estimationFunctions2.R")
source("Main Work/Code/simulationFunctions.R")

tt <- rep(Sys.time(), 2)
if(ncores > 1) requiredFunction <- c("Estimate.Loop", "Estimate.Loop2",
                                     "cor.matrix_to_norm.matrix", "triangle2vector","vector2triangle",
                                     "create_alpha_mat", "clean_sick", "vnorm", "compute_estimated_N","vector_var_matrix_calc_COR",
                                     "minusloglik", "bootstrapFunction")

ARMAdetails <- list(ARsick = c(0.4, -0.2), ARhealth = c(0.2, -0.1), 
                    MAsick = c(0.4), MAhealth = c(0.4))
sapply(ARMAdetails, checkInv)

Tlength <- 115
B <- 100
p <- 10
sampleDataB_ARMA <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                             percent_alpha = 0.4, range_alpha = c(0.6, 0.8), 
                             ARsick = ARMAdetails$ARsick , ARhealth = ARMAdetails$ARhealth,
                             MAsick = ARMAdetails$MAsick, MAhealth = ARMAdetails$MAhealth)

bootstrapFunction <- function(b){
  res_unspecified <- Estimate.Loop(Healthy_List = sampleDataB_ARMA$samples[[b]]$healthy,
                                   Sick_List = sampleDataB_ARMA$samples[[b]]$sick)
  
  res_specified <- Estimate.Loop2(theta0 = res_unspecified$theta,
                                  alpha0 = res_unspecified$alpha,
                                  healthy.data = sampleDataB_ARMA$samples[[b]]$healthy,
                                  sick.data = sampleDataB_ARMA$samples[[b]]$sick,
                                  T_thresh = 10^4, progress = FALSE)
  return(res_specified)
}

tt1 <- Sys.time()
if(ncores > 1){
  buildCL(ncores, c("dplyr", "matrixcalc"), requiredFunction)
  clusterExport(cl = cl, "sampleDataB_ARMA")
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

VarAlphaByHess <- ComputeFisher(simuldat[[1]], sampleDataB_ARMA$samples[[1]]$sick, "Hess") %>% solve
VarAlphaByGrad <- ComputeFisher(simuldat[[1]], sampleDataB_ARMA$samples[[1]]$sick, "Grad") %>% solve
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
  geom_histogram(bins = sqrt(B), col = "white", fill = "lightblue") + labs(title = "Loss of DF")
BiasRatio <- ggplot(data.frame(Bias = estN_all/Tlength - 1), aes(x = Bias)) +
  geom_histogram(bins = sqrt(B), col = "white", fill = "lightblue") + labs(title = "Loss of DF")

p <- 10
B <- 100
Tlength <- 115
MAlist <- c(-0.5, -0.25, 0, 0.1, 0.25, 0.4, 0.6)
lngth_MAlist <- length(MAlist)
seed <- sample(1:10000, 3)
sampleDataBT_ARMA <- list()

for(t in 1:lngth_MAlist){
  message(paste0("With MA = ", MAlist[t], "; ", t, "/", lngth_MAlist, " ", round(100*t/lngth_MAlist), "%"))
  sampleDataBT_ARMA[[t]] <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                                     MAsick = MAlist[t], MAhealth = MAlist[t],
                                     percent_alpha = 0.4, range_alpha = c(0.6, 0.8), seed = seed) %>%
    append(c("MA" = MAlist[t]), 0)
}

bootstrapFunction <- function(b, k){
  res_unspecified <- Estimate.Loop(Healthy_List = sampleDataBT_ARMA[[k]]$samples[[b]]$healthy,
                                   Sick_List = sampleDataBT_ARMA[[k]]$samples[[b]]$sick)
  
  res_specified <- Estimate.Loop2(theta0 = res_unspecified$theta,
                                  alpha0 = res_unspecified$alpha,
                                  healthy.data = sampleDataBT_ARMA[[k]]$samples[[b]]$healthy,
                                  sick.data = sampleDataBT_ARMA[[k]]$samples[[b]]$sick,
                                  T_thresh = 10^4, progress = FALSE)
  return(res_specified)
}

simuldatT <- list()

pb <- progress_bar$new(
  format = "Estimating models [:bar] :percent. Elapsed: :elapsed, ETA: :eta",
  total = lngth_MAlist, clear = FALSE, width= 90)

tt2 <- Sys.time()
if(ncores > 1){
  buildCL(ncores, c("dplyr", "matrixcalc"), requiredFunction)
  clusterExport(cl = cl, "sampleDataBT_ARMA")
  for(t in 1:lngth_MAlist){
    pb$tick()
    clusterExport(cl = cl, "t")
    simuldatT[[t]] <- parLapply(cl = cl, 1:B, bootstrapFunction, k = t)
  }
  terminateCL()
} else for(t in 1:lngth_MAlist){
  pb$tick()
  simuldatT[[t]] <- parLapply(cl = cl, 1:B, bootstrapFunction, k = t)
} 
tt2 <- Sys.time() - tt2
rm(pb)

alpha_simul <- array(dim = c(B, p, lngth_MAlist))
alpha_sdGrad <- matrix(0, nrow = lngth_MAlist, ncol = p)
alpha_sdHess <- matrix(0, nrow = lngth_MAlist, ncol = p)
alpha_sdComb <- matrix(0, nrow = lngth_MAlist, ncol = p)
emp_sds <- matrix(nrow = lngth_MAlist, ncol = p)
coeffs <- matrix(0, nrow = lngth_MAlist, ncol = 3)
estNT_all <- matrix(0, nrow = B, ncol = lngth_MAlist)

pb <- progress_bar$new(
  format = "Computing hessians [:bar] :percent. Elapsed: :elapsed, ETA: :eta",
  total = lngth_MAlist, clear = FALSE, width= 90)

for(t in 1:lngth_MAlist){
  for(b in 1:B) {
    alpha_simul[b,,t] <- simuldatT[[t]][[b]]$alpha
    estNT_all[b,t] <- simuldatT[[t]][[b]]$Est_N
  }
  pb$tick()
  tmpG <- ComputeFisher(simuldatT[[t]][[1]], sampleDataBT_ARMA[[t]]$samples[[1]]$sick, "Grad", silent = TRUE) %>% solve
  tmpH <- ComputeFisher(simuldatT[[t]][[1]], sampleDataBT_ARMA[[t]]$samples[[1]]$sick, "Hess", silent = TRUE) %>% solve
  tmpC <- tmpH %*% solve(tmpG) %*% tmpH
  
  alpha_sdGrad[t,] <- sqrt(diag(tmpG))
  alpha_sdHess[t,] <- sqrt(diag(tmpH))
  alpha_sdComb[t,] <- sqrt(diag(tmpC))
  
  emp_sds[t, ] <- apply(alpha_simul[,,t], 2, sd)
  coeffs[t, 1] <- lm(emp_sds[t,] ~ 0 + alpha_sdGrad[t,])$coef
  coeffs[t, 2] <- lm(emp_sds[t,] ~ 0 + alpha_sdHess[t,])$coef
  coeffs[t, 3] <- lm(emp_sds[t,] ~ 0 + alpha_sdComb[t,])$coef
}
rm(pb)

coeffs <- as.data.frame(coeffs)
colnames(coeffs) <- c("Grad", "Hess", "Combination")
coeffs$MA <- MAlist

CoefByDF <- gather(coeffs, key = Type, value = Coef, -MA) %>% ggplot(aes(x = MA, y = Coef, col = Type)) +
  geom_smooth(method = "lm", formula = y ~ (x + I(x^2)), se = FALSE) + geom_point()

ErrorByDF_Grad <- 
  inner_join(by = c("MAlist", "P"), cbind(MAlist, alpha_sdGrad) %>% as.data.frame() %>% gather(key = P, value = Value, -MAlist),
             cbind(MAlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -MAlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(MAlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Grad", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

ErrorByDF_Hess <- 
  inner_join(by = c("MAlist", "P"), cbind(MAlist, alpha_sdHess) %>% as.data.frame() %>% gather(key = P, value = Value, -MAlist),
             cbind(MAlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -MAlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(MAlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Hess", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

ErrorByDF_Combined <-
  inner_join(by = c("MAlist", "P"), cbind(MAlist, alpha_sdComb) %>% as.data.frame() %>% gather(key = P, value = Value, -MAlist),
             cbind(MAlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -MAlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(MAlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Combined", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

BiasDiffEstN <- as.data.frame(estNT_all - Tlength)
colnames(BiasDiffEstN) <- MAlist

EstNDiff <- gather(BiasDiffEstN, key = MA, value = Bias) %>%
  ggplot(aes(x = factor(MA, levels = MAlist, ordered = TRUE), y = Bias)) +
  geom_hline(yintercept = 0) + geom_boxplot(fill = "lightblue") + labs(x = "MA", y = "Absolute Loss of DF")

BiasRatioEstN <- as.data.frame(estNT_all / Tlength - 1)
colnames(BiasRatioEstN) <- MAlist

EstNRatio <- gather(BiasRatioEstN, key = MA, value = Bias) %>%
  ggplot(aes(x = factor(MA, levels = MAlist, ordered = TRUE), y = Bias)) +
  geom_hline(yintercept = 0) + geom_boxplot(fill = "lightblue") + labs(x = "MA", y = "Relative Loss of DF")


link2 <- gsub(":", "-", paste0("Main Work/Data/Enviroments/", "fullRunYesARMA ", Sys.time(), ".RData") )

save.image(file = link2)


ErrorByDF_Combined
ErrorByDF_Grad
ErrorByDF_Hess

SDErrorByFunction
CoefByDF

BiasDiff
BiasRatio

EstNDiff
EstNRatio
