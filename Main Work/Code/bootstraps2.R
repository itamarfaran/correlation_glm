source("Main Work/Code/generalFunctions.R")
source("Main Work/code/estimationFunctions2.R")
source("Main Work/code/simulationFunctions.R")

tt <- rep(Sys.time(), 2)
ncores <- detectCores() - 1 #### Enter here!
if(ncores > 1) requiredFunction <- c("Estimate.Loop", "Estimate.Loop2",
                                     "cor.matrix_to_norm.matrix", "triangle_to_vector","vector_to_triangle",
                                     "create_alpha_mat", "clean_sick", "vnorm", "compute_estimated_N","vector_var_matrix_calc_COR",
                                     "minusloglik", "bootstrapFunction")
Tlength <- 115
B <- 100
p <- 10
sampleData <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                            percent_alpha = 0.4, range_alpha = c(0.6, 0.8))


bootstrapFunction <- function(b){
  res_unspecified <- Estimate.Loop(Healthy_List = sampleData$samples[[b]]$healthy,
                                   Sick_List = sampleData$samples[[b]]$sick)
  
  res_specified <- Estimate.Loop2(theta0 = res_unspecified$theta,
                                  alpha0 = res_unspecified$alpha,
                                  healthy.data = sampleData$samples[[b]]$healthy,
                                  sick.data = sampleData$samples[[b]]$sick,
                                  T_thresh = 10^4, progress = FALSE)
  return(res_specified)
}

tt[1] <- Sys.time()
if(ncores > 1){
  buildcl(ncores, c("dplyr", "matrixcalc"), requiredFunction)
  clusterExport(cl = cl, "sampleData")
  simuldat <- parLapply(cl = cl, 1:B, bootstrapFunction)
  stopCluster(cl = cl)
} else {
  simuldat <- lapply(1:B, bootstrapFunction)
}
tt[1] <- Sys.time() - tt[1]
alpha_simul <- matrix(nrow = B, ncol = p)
estN_all <- numeric(B)

for(b in 1:B){
  alpha_simul[b,] <- simuldat[[b]]$alpha
  estN_all[b] <- simuldat[[b]]$Est_N
} 
VarAlphaByHess <- ComputeFisher(simuldat[[1]], sampleData$samples[[1]]$sick, "Hess") %>% solve
VarAlphaByGrad <- ComputeFisher(simuldat[[1]], sampleData$samples[[1]]$sick, "Grad") %>% solve
VarAlphaCombined <- VarAlphaByHess %*% solve(VarAlphaByGrad) %*% VarAlphaByHess

Emp_vs_Theo <- data.frame(TheoreticHess = sqrt(diag(VarAlphaByHess)),
                          TheoreticGrad = sqrt(diag(VarAlphaByGrad)),
                          TheoreticCombined = sqrt(diag(VarAlphaByHess %*% solve(VarAlphaByGrad) %*% VarAlphaByHess)),
                          Empiric = sapply(1:p, function(i) sd(alpha_simul[,i]))) %>%
  mutate(QuotentHess = TheoreticHess/Empiric,
         QuotentGrad = TheoreticGrad/Empiric,
         QuotentCombined = TheoreticCombined/Empiric,
         QuotentBetween = TheoreticHess/TheoreticGrad)
  

Tlength
estN_all
estN_all/Tlength

link2 <- paste0("Main Work/Data/Enviroments/enviroment ",
                "p", p, " ", "B", B, " ", gsub(":", "", Sys.time()), ".RData")
save.image(link2)
rm(link2)



p <- 6
B <- 100
Tlist <- c(10, 30, 50, 70, 100, 120, 150, 170, 200, 250, 300, 400, 700, 1000, 1500, 2000, 3000, 4000)
lngth_Tlist <- length(Tlist)
seed <- sample(1:10000, 3)
sampleDataT <- list()

for(t in 1:lngth_Tlist){
  sampleDataT[[t]] <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlist[t],
                                    percent_alpha = 0.4, range_alpha = c(0.6, 0.8), seed = seed) %>%
    append(c("Tlength" = Tlist[t]), 0)
}


bootstrapFunction <- function(b, k){
  res_unspecified <- Estimate.Loop(Healthy_List = sampleDataT[[k]]$samples[[b]]$healthy,
                                   Sick_List = sampleDataT[[k]]$samples[[b]]$sick)
  
  res_specified <- Estimate.Loop2(theta0 = res_unspecified$theta,
                                  alpha0 = res_unspecified$alpha,
                                  healthy.data = sampleDataT[[k]]$samples[[b]]$healthy,
                                  sick.data = sampleDataT[[k]]$samples[[b]]$sick,
                                  T_thresh = 10^4, progress = FALSE)
  return(res_specified)
}
simuldatT <- list()

tt[2] <- Sys.time()
if(ncores > 1){
  buildcl(ncores, c("dplyr", "matrixcalc"), requiredFunction)
  clusterExport(cl = cl, "sampleDataT")
  for(t in 1:lngth_Tlist){
    cat(paste0(Tlist[t], " (", round(100*t/lngth_Tlist), "%); "))
    clusterExport(cl = cl, "t")
    simuldatT[[t]] <- parLapply(cl = cl, 1:B, bootstrapFunction, k = t)
  }
  stopCluster(cl = cl)
} else for(t in 1:lngth_Tlist){
  cat(paste0(Tlist[t], " (", round(100*t/lngth_Tlist), "%); "))
  simuldatT[[t]] <- parLapply(cl = cl, 1:B, bootstrapFunction, k = t)
} 
tt[2] <- Sys.time() - tt[2]

alpha_simul <- array(dim = c(B, p, lngth_Tlist))
alpha_sdGrad <- matrix(0, nrow = lngth_Tlist, ncol = p)
alpha_sdHess <- matrix(0, nrow = lngth_Tlist, ncol = p)
alpha_sdComb <- matrix(0, nrow = lngth_Tlist, ncol = p)
emp_sds <- matrix(nrow = lngth_Tlist, ncol = p)
coeffs <- matrix(0, nrow = lngth_Tlist, ncol = 3)
estNT_all <- matrix(0, nrow = B, ncol = lngth_Tlist)

for(t in 1:lngth_Tlist){
  for(b in 1:B) {
    alpha_simul[b,,t] <- simuldatT[[t]][[b]]$alpha
    estNT_all[b,t] <- simuldatT[[t]][[b]]$Est_N
  }
  tmpG <- ComputeFisher(simuldatT[[t]][[1]], sampleDataT[[t]]$samples[[1]]$sick, "Grad") %>% solve
  tmpH <- ComputeFisher(simuldatT[[t]][[1]], sampleDataT[[t]]$samples[[1]]$sick, "Hess") %>% solve
  tmpC <- tmpH %*% solve(tmpG) %*% tmpH
  
  alpha_sdGrad[t,] <- sqrt(diag(tmpG))
  alpha_sdHess[t,] <- sqrt(diag(tmpH))
  alpha_sdComb[t,] <- sqrt(diag(tmpC))
  
  emp_sds[t, ] <- apply(alpha_simul[,,t], 2, sd)
  coeffs[t, 1] <- lm(emp_sds[t,] ~ 0 + alpha_sdGrad[t,])$coef
  coeffs[t, 2] <- lm(emp_sds[t,] ~ 0 + alpha_sdHess[t,])$coef
  coeffs[t, 3] <- lm(emp_sds[t,] ~ 0 + alpha_sdComb[t,])$coef
}


sdsDims <- dim(emp_sds)
forLm <- matrix(0, nrow = prod(sdsDims), ncol = 4)
for(i in 1:sdsDims[1]){
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 1] <- Tlist[i]
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 2] <- alpha_sd_est[i,]
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 3] <- alpha_sd_est2[i,]
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 4] <- emp_sds[i,]
}

forLm <- as.data.frame(forLm)
colnames(forLm) <- c("DF", "EstimatedHess", "EstimatedGrad", "Empiric")
forLm <- mutate(forLm, QuotentHess = EstimatedHess/Empiric,
                QuotentGrad = EstimatedGrad/Empiric,
                QuotentBetween = EstimatedGrad/EstimatedHess)

summary(aov(QuotentGrad ~ factor(DF), data = forLm))
summary(aov(QuotentBetween ~ factor(DF), data = forLm))
summary(lm(QuotentBetween ~ log(DF), data = forLm))

summary(lm(QuotentGrad ~ 0 + log(DF)*QuotentHess, data = forLm))


ggplot(forLm, aes(x = factor(DF), col = factor(DF), y = QuotentGrad)) + geom_point()
ggplot(forLm, aes(x = factor(DF), col = factor(DF), y = QuotentHess)) + geom_point()
ggplot(forLm, aes(x = factor(DF), col = factor(DF), y = QuotentBetween)) + geom_point()

tmp <- numeric()
for(i in 1:lngth_Tlist) tmp <- c(tmp, rep(Tlist[i], p))
forplt <- data.frame(DF = factor(tmp), Empiric = as.vector(t(emp_sds)), Estimate = as.vector(t(alpha_sd_est2)))
tmp <- rbind(eff_n, round(eff_n/Tlist, 3), coeffs)
colnames(tmp) <- Tlist
row.names(tmp) <- c("Est_n", "Ratio", "Coeffs")


ggplot(forplt, aes(x = Estimate, y = Empiric)) + geom_point(aes(col = DF)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, size = 0.8, linetype = 2, col = "darkgrey") +
  geom_smooth(method = "lm", formula = y ~ 0 + x + exp(x) ,se = FALSE, linetype = 4)
  # geom_abline(slope = 2, intercept = 0, size = 0.8, col = "darkred", linetype = 2)


data.frame(DF = Tlist, Coefficient = coeffs) %>% ggplot(aes(x = DF, y = Coefficient)) + geom_point() +
  geom_hline(yintercept = 1, col = "darkgrey") + geom_vline(xintercept = 0) + 
  geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE)

#View(t(tmp))
x <- 1/Tlist
summary(lm(coeffs ~ x))

df_error <- tmp[2,]-1
hist(df_error)
summary(df_error)
plot(Tlist, df_error)
summary(lm(df_error ~ 0 + Tlist))

link2 <- paste0("Main Work/Data/Enviroments/enviroment robustness_check p_",
                p, " T_", lngth_Tlist, " ", gsub(":", "", Sys.time()), ".RData")
save.image(link2)
rm(link2)

