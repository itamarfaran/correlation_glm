source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

tt.all <- Sys.time()
# profvis({

Tlength <- 115
ARMAdetails <- list(ARsick = 0.3, MAsick = NULL,
                    ARhealth = 0.3, MAhealth = NULL)
sapply(ARMAdetails, checkInv)
sampleData <- createSamples(nH = 57, nS = 42, p = 22, Tlength = Tlength,
                            percent_alpha = 0.3, range_alpha = c(0.65, 0.95),
                            ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
                            ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
                            ncores = ncores)

#Are all matrices positive definite?
all(abind(sampleData$healthy, sampleData$sick, along = 3) %>%
      apply(3, is.positive.definite))


Pelet_IID <- Estimate.Loop(sampleData$healthy, sampleData$sick, MaxLoop = 100)

tt.est <- Sys.time()
Pelet_Cov <- Estimate.Loop2(theta0 = Pelet_IID$theta, alpha0 = Pelet_IID$alpha,
                            healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                            T_thresh = Tlength, method = "Nelder-Mead", updateU = 1, progress = T)
tt.est <- Sys.time() - tt.est

tt.hess <- Sys.time()
fisherMatrHess <- ComputeFisher(Pelet_Cov, sampleData$sick, "Hess")  %>% regularizeMatrix()
tt.hess <- Sys.time() - tt.hess

tt.grad <- Sys.time()
fisherMatrGrad <- ComputeFisher(Pelet_Cov, sampleData$sick, "Grad", ncores = ncores)  %>% regularizeMatrix()
tt.grad <- Sys.time() - tt.grad

fisherMatrComb <- fisherMatrHess %*% solve(fisherMatrGrad) %*% fisherMatrHess

HypTestResHess <- build_hyp.test(Pelet_Cov, fisherMatrHess, sampleData$alpha, Real = sampleData$alpha)
HypTestResGrad <- build_hyp.test(Pelet_Cov, fisherMatrGrad, sampleData$alpha, Real = sampleData$alpha)
HypTestResComb <- build_hyp.test(Pelet_Cov, fisherMatrComb, sampleData$alpha, Real = sampleData$alpha)
gc()

Pelet_Cov$returns
Pelet_Cov$convergence
c("Est_DF" = Pelet_Cov$Est_N, "Real_DF" = Tlength)
c("Test" = HypTestResHess$Test, "Sig Level" = HypTestResHess$Significance, "FWER Method" = HypTestResHess$MH_method)

tt.hyp <- Sys.time()
HypTestResHess$Results[order(HypTestResHess$Results$Real),]
tt.hyp <- Sys.time() - tt.hyp
HypTestResGrad$Results[order(HypTestResGrad$Results$Real),]
HypTestResComb$Results[order(HypTestResComb$Results$Real),]

wilksTest(Pelet_Cov, sampleData$healthy, sampleData$sick)

tt.all <- Sys.time() - tt.all

# })