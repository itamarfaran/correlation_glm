source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$Exponent

tt <- rep(Sys.time(), 2)
ARMAdetails <- list(ARsick = c(0.4, -0.2), ARhealth = c(0.2, -0.1), 
                    MAsick = c(0.4), MAhealth = c(0.4))
sapply(ARMAdetails, checkInv)

Tlength <- 115
B <- 120
p <- 32
sampleDataB_ARMA <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                             percent_alpha = 0.4, range_alpha = c(0.6, 0.8), 
                             ARsick = ARMAdetails$ARsick , ARhealth = ARMAdetails$ARhealth,
                             MAsick = ARMAdetails$MAsick, MAhealth = ARMAdetails$MAhealth, ncores = ncores)

bootstrapFunction <- function(b) estimateAlpha(healthy.data = sampleDataB_ARMA$samples[[b]]$healthy,
                                               sick.data = sampleDataB_ARMA$samples[[b]]$sick,
                                               T_thresh = 10^4, updateU = 1, progress = F)


tt1 <- Sys.time()
simuldat <- mclapply(1:B, bootstrapFunction, mc.cores = ncores)
tt1 <- Sys.time() - tt1

alpha_simul <- matrix(nrow = B, ncol = p)
estN_all <- numeric(B)

for(b in 1:B){
  alpha_simul[b,] <- simuldat[[b]]$alpha
  estN_all[b] <- simuldat[[b]]$Est_N
}

VarAlphaByHess <- ComputeFisher(simuldat[[1]], sampleDataB_ARMA$samples[[1]]$sick, "Hess") %>% solve
VarAlphaByGrad <- ComputeFisher(simuldat[[1]], sampleDataB_ARMA$samples[[1]]$sick, "Grad", ncores = ncores) %>% solve
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

p <- 32
B <- 100
Tlength <- 100
ARlist <- c(-0.5, -0.25, 0, 0.1, 0.25, 0.4, 0.6)
lngth_ARlist <- length(ARlist)
seed <- sample(1:10000, 3)
sampleDataBT_ARMA <- list()

for(t in 1:lngth_ARlist){
  # cat(paste0("With AR = ", ARlist[t], "; ", t, "/", lngth_ARlist, " ", round(100*t/lngth_ARlist), "%; "))
  sampleDataBT_ARMA[[t]] <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                                     ARsick = ARlist[t], ARhealth = ARlist[t],
                                     percent_alpha = 0.4, range_alpha = c(0.6, 0.8), seed = seed, ncores = ncores) %>%
    append(c("AR" = ARlist[t]), 0)
}


bootstrapFunction <- function(b, k) estimateAlpha(healthy.data = sampleDataBT_ARMA[[k]]$samples[[b]]$healthy,
                                               sick.data = sampleDataBT_ARMA[[k]]$samples[[b]]$sick,
                                               T_thresh = 10^4, updateU = 1, progress = F)

simuldatT <- list()

pb <- progress_bar$new(
  format = "Estimating models [:bar] :percent. Elapsed: :elapsed, ETA: :eta",
  total = lngth_ARlist, clear = FALSE, width= 90)

tt2 <- Sys.time()
for(t in 1:lngth_ARlist){
  pb$tick()
  simuldatT[[t]] <- mclapply(1:B, bootstrapFunction, k = t, mc.cores = ncores)
}
tt2 <- Sys.time() - tt2
rm(pb)

alpha_simul <- array(dim = c(B, p, lngth_ARlist))
alpha_sdGrad <- matrix(0, nrow = lngth_ARlist, ncol = p)
alpha_sdHess <- matrix(0, nrow = lngth_ARlist, ncol = p)
alpha_sdComb <- matrix(0, nrow = lngth_ARlist, ncol = p)
emp_sds <- matrix(nrow = lngth_ARlist, ncol = p)
coeffs <- matrix(0, nrow = lngth_ARlist, ncol = 3)
estNT_all <- matrix(0, nrow = B, ncol = lngth_ARlist)

pb <- progress_bar$new(
  format = "Computing Fisher information [:bar] :percent. Elapsed: :elapsed, ETA: :eta",
  total = lngth_ARlist, clear = FALSE, width= 90)

for(t in 1:lngth_ARlist){
  for(b in 1:B) {
    alpha_simul[b,,t] <- simuldatT[[t]][[b]]$alpha
    estNT_all[b,t] <- simuldatT[[t]][[b]]$Est_N
  }
  pb$tick()
  tmpG <- ComputeFisher(simuldatT[[t]][[1]], sampleDataBT_ARMA[[t]]$samples[[1]]$sick, "Grad", silent = TRUE, ncores = ncores) %>% solve
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
coeffs$AR <- ARlist

CoefByDF <- gather(coeffs, key = Type, value = Coef, -AR) %>% ggplot(aes(x = AR, y = Coef, col = Type)) +
  geom_smooth(method = "lm", formula = y ~ (x + I(x^2)), se = FALSE) + geom_point()

ErrorByDF_Grad <- 
  inner_join(by = c("ARlist", "P"), cbind(ARlist, alpha_sdGrad) %>% as.data.frame() %>% gather(key = P, value = Value, -ARlist),
             cbind(ARlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -ARlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(ARlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Grad", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

ErrorByDF_Hess <- 
  inner_join(by = c("ARlist", "P"), cbind(ARlist, alpha_sdHess) %>% as.data.frame() %>% gather(key = P, value = Value, -ARlist),
             cbind(ARlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -ARlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(ARlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Hess", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

ErrorByDF_Combined <-
  inner_join(by = c("ARlist", "P"), cbind(ARlist, alpha_sdComb) %>% as.data.frame() %>% gather(key = P, value = Value, -ARlist),
             cbind(ARlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -ARlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(ARlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Combined", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

BiasDiffEstN <- as.data.frame(estNT_all - Tlength)
colnames(BiasDiffEstN) <- ARlist

EstNDiff <- gather(BiasDiffEstN, key = AR, value = Bias) %>%
  ggplot(aes(x = factor(AR, levels = ARlist, ordered = TRUE), y = Bias)) +
  geom_hline(yintercept = 0) + geom_boxplot(fill = "lightblue") + labs(x = "AR", y = "Absolute Loss of DF")

BiasRatioEstN <- as.data.frame(estNT_all / Tlength - 1)
colnames(BiasRatioEstN) <- ARlist

EstNRatio <- gather(BiasRatioEstN, key = AR, value = Bias) %>%
  ggplot(aes(x = factor(AR, levels = ARlist, ordered = TRUE), y = Bias)) +
  geom_hline(yintercept = 0) + geom_boxplot(fill = "lightblue") + labs(x = "AR", y = "Relative Loss of DF")


link2 <- gsub(":", "-", paste0("main_work/Data/Enviroments/", "fullRunYesARMA ", Sys.time(), ".RData") )

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
