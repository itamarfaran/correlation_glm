source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$Exponent

Tlength <- 115
B <- 120
p <- 32
sampleDataB <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                             percent_alpha = 0.4, range_alpha = c(0.6, 0.8), ncores = ncores)


bootstrapFunction <- function(b) estimateAlpha(healthy.data = sampleDataB$samples[[b]]$healthy,
                                               sick.data = sampleDataB$samples[[b]]$sick,
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
VarAlphaByHess <- ComputeFisher(simuldat[[1]], sampleDataB$samples[[1]]$sick, "Hess") %>% solve
VarAlphaByGrad <- ComputeFisher(simuldat[[1]], sampleDataB$samples[[1]]$sick, "Grad", ncores = ncores) %>% solve
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

p <- 32
B <- 80
Tlist <- c(10, 30, 50, 70, 100, 120, 150, 170, 200, 250, 300, 400, 700, 1000, 1500, 2000, 3000, 4000)
lngth_Tlist <- length(Tlist)
seed <- sample(1:10000, 3)
sampleDataBT <- list()

for(t in 1:lngth_Tlist){
  sampleDataBT[[t]] <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlist[t],
                                     percent_alpha = 0.4, range_alpha = c(0.6, 0.8), seed = seed, ncores = ncores) %>%
    append(c("Tlength" = Tlist[t]), 0)
}


bootstrapFunction <- function(b, k) estimateAlpha(healthy.data = sampleDataBT[[k]]$samples[[b]]$healthy,
                                                  sick.data = sampleDataBT[[k]]$samples[[b]]$sick,
                                                  T_thresh = 10^4, updateU = 1, progress = F)

simuldatT <- list()

pb <- progress_bar$new(
  format = "Estimating models [:bar] :percent. Elapsed: :elapsed, ETA: :eta",
  total = lngth_Tlist, clear = FALSE, width= 90)

tt2 <- Sys.time()
for(t in 1:lngth_Tlist){
  simuldatT[[t]] <- mclapply(1:B, bootstrapFunction, k = t, mc.cores = ncores)
  pb$tick()
}
tt2 <- Sys.time() - tt2

alpha_simul <- array(dim = c(B, p, lngth_Tlist))
alpha_sdGrad <- matrix(0, nrow = lngth_Tlist, ncol = p)
alpha_sdHess <- matrix(0, nrow = lngth_Tlist, ncol = p)
alpha_sdComb <- matrix(0, nrow = lngth_Tlist, ncol = p)
emp_sds <- matrix(nrow = lngth_Tlist, ncol = p)
coeffs <- matrix(0, nrow = lngth_Tlist, ncol = 3)
estNT_all <- matrix(0, nrow = B, ncol = lngth_Tlist)

pb <- progress_bar$new(
  format = "Computing Fisher information [:bar] :percent. Elapsed: :elapsed, ETA: :eta",
  total = lngth_Tlist, clear = FALSE, width= 90)

for(t in 1:lngth_Tlist){
  for(b in 1:B) {
    alpha_simul[b,,t] <- simuldatT[[t]][[b]]$alpha
    estNT_all[b,t] <- simuldatT[[t]][[b]]$Est_N
  }
  pb$tick()
  tmpG <- ComputeFisher(simuldatT[[t]][[1]], sampleDataBT[[t]]$samples[[1]]$sick, "Grad", silent = TRUE, ncores = ncores) %>% solve
  tmpH <- ComputeFisher(simuldatT[[t]][[1]], sampleDataBT[[t]]$samples[[1]]$sick, "Hess", silent = TRUE) %>% solve
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
coeffs$Tlength <- Tlist

CoefByDF <- gather(coeffs, key = Type, value = Coef, -Tlength) %>% ggplot(aes(x = Tlength, y = Coef, col = Type)) +
  geom_smooth(se = FALSE) + geom_point()

ErrorByDF_Grad <- 
  inner_join(by = c("Tlist", "P"), cbind(Tlist, alpha_sdGrad) %>% as.data.frame() %>% gather(key = P, value = Value, -Tlist),
             cbind(Tlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -Tlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(Tlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Grad", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

ErrorByDF_Hess <- 
  inner_join(by = c("Tlist", "P"), cbind(Tlist, alpha_sdHess) %>% as.data.frame() %>% gather(key = P, value = Value, -Tlist),
             cbind(Tlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -Tlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(Tlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Hess", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

ErrorByDF_Combined <-
  inner_join(by = c("Tlist", "P"), cbind(Tlist, alpha_sdComb) %>% as.data.frame() %>% gather(key = P, value = Value, -Tlist),
             cbind(Tlist, emp_sds) %>% as.data.frame() %>% gather(key = P, value = Value, -Tlist) ) %>%
  ggplot(aes(x = Value.x, y = Value.y, col = factor(Tlist))) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_abline(slope = 1, intercept = 0, col = "blue", linetype = 2, size = 1) +
  geom_point() + labs(x = "Theoritcal by Combined", y = "Empiric") + xlim(0, 0.15) + ylim(0, 0.15)

BiasDiffEstN <- as.data.frame(estNT_all - rep(1, B) %*% t(Tlist))
colnames(BiasDiffEstN) <- Tlist

EstNDiff <- gather(BiasDiffEstN, key = DF, value = Bias) %>%
  ggplot(aes(x = factor(DF, levels = Tlist, ordered = TRUE), y = Bias)) +
  geom_hline(yintercept = 0) + geom_boxplot(fill = "lightblue") + labs(x = "DF", y = "Absolute Bias")

BiasRatioEstN <- as.data.frame(estNT_all / rep(1, B) %*% t(Tlist) - 1)
colnames(BiasRatioEstN) <- Tlist

EstNRatio <- gather(BiasRatioEstN, key = DF, value = Bias) %>%
  ggplot(aes(x = factor(DF, levels = Tlist, ordered = TRUE), y = Bias)) +
  geom_hline(yintercept = 0) + geom_boxplot(fill = "lightblue") + labs(x = "DF", y = "Relative Bias")


link2 <- gsub(":", "-", paste0("main_work/Data/Enviroments/", "fullRunNoARMA ", Sys.time(), ".RData") )

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
