source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

B = 50

linkFun <- linkFunctions$multiplicative_identity

ARMAdetails <- list(ARsick = c(0.4, -0.2), ARhealth = c(0.2, -0.1), 
                    MAsick = c(0.4), MAhealth = c(0.4))
sapply(ARMAdetails, checkInv)


compare_fdr_and_pwr <- function(B = 40, nH = 19, nS = 12, p = 10, percent_alpha = 0.3,
                                range_alpha = c(0.8, 0.95), Tlength = 115,
                                linkFun = linkFunctions$multiplicative_identity,
                                ARMAdetails = list(ARsick = 0, ARhealth = 0, 
                                                   MAsick = 0, MAhealth = 0),
                                ncores = 1){
  sampleDataB_ARMA <-
    createSamples(B = B, nH = nH, nS = nS, p = p, Tlength = Tlength,
                  percent_alpha = percent_alpha, range_alpha = range_alpha,
                  ARsick = ARMAdetails$ARsick , ARhealth = ARMAdetails$ARhealth,
                  MAsick = ARMAdetails$MAsick, MAhealth = ARMAdetails$MAhealth,
                  ncores = ncores)
  
  estimation_results <- mclapply(1:B, function(b){
    estimateAlpha(healthy.data = sampleDataB_ARMA$samples[[b]]$healthy,
                  sick.data = sampleDataB_ARMA$samples[[b]]$sick,
                  linkFun = linkFun, T_thresh = 10^4, updateU = 1, progress = F)}, mc.cores = ncores)
  estimated_variances <- simplify2array(mclapply(1:B, function(b){
    FisherByHess <- ComputeFisher(estimation_results[[b]], sampleDataB_ARMA$samples[[b]]$sick, "Hess", linkFun = linkFun)
    FisherByGrad <- ComputeFisher(estimation_results[[b]], sampleDataB_ARMA$samples[[b]]$sick, "Grad", linkFun = linkFun, ncores = 1)
    SolveFisherByHess <- solve(FisherByHess)
    VarAlphaCombined <- SolveFisherByHess %*% FisherByGrad %*% SolveFisherByHess
    return(VarAlphaCombined)
  }, mc.cores = ncores))
  alphas <- do.call(rbind, transpose(estimation_results)$alpha)
  alphas_sd <- t(sapply(1:B, function(b) sqrt(diag(estimated_variances[,,b]))))
  alphas_pval <- 2*pnorm(abs((alphas - linkFun$NULL_VAL)/alphas_sd), lower.tail = F)
  ttest_res <- do.call(rbind, mclapply(1:B, function(b) { multipleComparison(
    sampleDataB_ARMA$samples[[b]]$healthy, sampleDataB_ARMA$samples[[b]]$sick,
    p.adjust.method = "none", test = "both")}, mc.cores = ncores))
  
  which_alphas_rejected <- t(apply(alphas_pval, 1, p.adjust, method = "holm") < 0.05)
  ttest_res_matrices <- sapply(1:B, function(b) vector2triangle(ttest_res[b,]), simplify = "array")
  twostep_hyptest <- do.call(rbind, mclapply(1:B, function(b){
    tt_alpha = which_alphas_rejected[b,]
    tt_pvals = ttest_res_matrices[,,b]
    
    tt_pvals[,!tt_alpha] <- 1
    if(sum(tt_alpha) > 0) tt_pvals[,tt_alpha] <- apply(tt_pvals[,tt_alpha], 2, p.adjust, method = "BH")
    triangle2vector(tt_pvals)
  }))
  classic_hyptest <- t(apply(ttest_res, 1, p.adjust, method = "BH"))
  
  real_correlation_changes <- triangle2vector(linkFun$FUN(rep(1, p*(p-1)/2), sampleDataB_ARMA$alpha, 1))
  null_alphas <- real_correlation_changes == 1
  
  two_step_FDR <- mean(rowMeans(twostep_hyptest[, which(null_alphas)] < 0.05))
  classic_FDR <- mean(rowMeans(classic_hyptest[, which(null_alphas)] < 0.05))
  
  two_step_power <- mean(rowMeans(twostep_hyptest[, which(!null_alphas)] < 0.05))
  classic_power <- mean(rowMeans(classic_hyptest[, null_alphas] < 0.05))
  
  output <- c(two_step_FDR, classic_FDR, two_step_power, classic_power)
  names(output) <- c("two_step_FDR", "classic_FDR", "two_step_power", "classic_power")
  
  link <- gsub(":", "-", paste0("main_work/Data/Enviroments/", "pwr_comp_innerenv_", Sys.time(), ".RData"))
  save.image(file = link)
  return(output)
}


experiment_1 <- # Fixed p, increasing effect size
  sapply(list(c(0.8, 0.85), c(0.85, 0.9), c(0.9, 0.95)),
         function(interval) compare_fdr_and_pwr(
           B = B, nH = 40, nS = 40, p = 20, percent_alpha = 0.3,
           range_alpha = interval, Tlength = 115,
           linkFun = linkFunctions$multiplicative_identity,
           ARMAdetails = ARMAdetails, ncores = ncores))

link2 <- gsub(":", "-", paste0("main_work/Data/Enviroments/", "power_comp ", Sys.time(), ".RData") )
save.image(file = link2)

experiment_2 <- # Fixed alpha, increasing dimension
  sapply(5*(2:8), function(p) compare_fdr_and_pwr(
    B = B, nH = 60, nS = 60, p = p, percent_alpha = 0.3,
    range_alpha = c(0.8, 0.8), Tlength = 115,
    linkFun = linkFunctions$multiplicative_identity,
    ARMAdetails = ARMAdetails, ncores = ncores))

file.remove(link2)
link2 <- gsub(":", "-", paste0("main_work/Data/Enviroments/", "power_comp ", Sys.time(), ".RData") )
save.image(file = link2)

experiment_3 <- # Fixed alpha, increasing dimension
  sapply(5*(3:8), function(p) compare_fdr_and_pwr(
    B = B, nH = 60, nS = 60, p = p, percent_alpha = 0.3,
    range_alpha = c(0.9, 0.9), Tlength = 115,
    linkFun = linkFunctions$multiplicative_identity,
    ARMAdetails = ARMAdetails, ncores = ncores))

file.remove(link2)
link2 <- gsub(":", "-", paste0("main_work/Data/Enviroments/", "power_comp ", Sys.time(), ".RData") )
save.image(file = link2)
