source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

tt <- Sys.time()
# profvis({

Tlength <- 115
ARMAdetails <- list(ARsick = 0.3, MAsick = NULL,
                    ARhealth = 0.3, MAhealth = NULL)
sapply(ARMAdetails, checkInv)
sampleData <- createSamples(nH = 57, nS = 42, p = 10, Tlength = Tlength,
                                percent_alpha = 0.3, range_alpha = c(0.65, 0.95),
                            ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
                            ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth)

#Are all matrices positive definite?
all(abind(sampleData$healthy, sampleData$sick, along = 3) %>%
      apply(3, is.positive.definite))


Pelet_IID <- Estimate.Loop(sampleData$healthy, sampleData$sick, MaxLoop = 100)

# emp <- sampleData$healthy %>% cor.matrix_to_norm.matrix() %>% cov() %>% triangle2vector(diag = TRUE)
# theo <- sampleData$healthy %>% calculate_mean_matrix() %>% vector_var_matrix_calc_COR() %>%
#   triangle2vector(diag = TRUE)
# 
# data.frame(Theoretical = theo, Empirical = emp) %>% ggplot(aes(x = Theoretical, y = Empirical)) + 
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   geom_point(col = "blue", alpha = 0.5) + geom_smooth(method = "lm")

Pelet_Cov <- Estimate.Loop2(theta0 = Pelet_IID$theta, alpha0 = Pelet_IID$alpha,
                            healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                            T_thresh = Tlength, method = "Nelder-Mead")

fisherMatrHess <- ComputeFisher(Pelet_Cov, sampleData$sick, "Hess")  %>% regularizeMatrix()
fisherMatrGrad <- ComputeFisher(Pelet_Cov, sampleData$sick, "Grad")  %>% regularizeMatrix()
fisherMatrComb <- fisherMatrHess %*% solve(fisherMatrGrad) %*% fisherMatrHess

HypTestResHess <- build_hyp.test(Pelet_Cov, fisherMatrHess, sampleData$alpha, Real = sampleData$alpha)
HypTestResGrad <- build_hyp.test(Pelet_Cov, fisherMatrGrad, sampleData$alpha, Real = sampleData$alpha)
HypTestResComb <- build_hyp.test(Pelet_Cov, fisherMatrComb, sampleData$alpha, Real = sampleData$alpha)
gc()

Pelet_Cov$returns
Pelet_Cov$convergence
c("Est_DF" = Pelet_Cov$Est_N, "Real_DF" = Tlength)
c("Test" = HypTestResHess$Test, "Sig Level" = HypTestResHess$Significance, "FWER Method" = HypTestResHess$MH_method)

HypTestResHess$Results[order(HypTestResHess$Results$Real),]
HypTestResGrad$Results[order(HypTestResGrad$Results$Real),]
HypTestResComb$Results[order(HypTestResComb$Results$Real),]

#for(i in 2:Pelet_Cov$returns) print(Pelet_Cov$Log_Optim[[i]]$counts)

wilksTest(Pelet_Cov, sampleData$healthy, sampleData$sick)

tt <- Sys.time() - tt

# })
tt