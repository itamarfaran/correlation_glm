source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

linkFun <- linkFunctions$Identity

tt.all <- Sys.time()
# profvis({

Tlength <- 115
ARMAdetails <- list(ARsick = 0.3, MAsick = NULL,
                    ARhealth = 0.3, MAhealth = NULL)
sapply(ARMAdetails, checkInv)
sampleData <- createSamples(nH = 87, nS = 63, p = 12, Tlength = Tlength, dim_alpha = 1,
                            percent_alpha = 0.3, range_alpha = c(0.65, 0.95),
                            ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
                            ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth,
                            ncores = ncores)

#Are all matrices positive definite?
all(abind(sampleData$healthy, sampleData$sick, along = 3) %>%
      apply(3, is.positive.definite))

tt.est <- Sys.time()
Pelet_Cov <- estimateAlpha(healthy.data = sampleData$healthy, sick.data = sampleData$sick, dim_alpha = 1,
                            T_thresh = Tlength, updateU = 1, progress = T, linkFun = linkFun)
tt.est <- Sys.time() - tt.est


tt.hess <- Sys.time()
fisherMatrHess <- ComputeFisher(Pelet_Cov, sampleData$sick, "Hess", linkFun = linkFun, dim_alpha = 1)  %>% regularizeMatrix()
tt.hess <- Sys.time() - tt.hess

tt.grad <- Sys.time()
fisherMatrGrad <- ComputeFisher(Pelet_Cov, sampleData$sick, "Grad", linkFun = linkFun, ncores = ncores, dim_alpha = 1)  %>% regularizeMatrix()
tt.grad <- Sys.time() - tt.grad

fisherMatrSandwich <- fisherMatrHess %*% forceSolve(fisherMatrGrad) %*% fisherMatrHess

HypTestResHess <- build_hyp.test(Pelet_Cov, fisherMatrHess, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResGrad <- build_hyp.test(Pelet_Cov, fisherMatrGrad, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResComb <- build_hyp.test(Pelet_Cov, fisherMatrSandwich, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
gc()

Pelet_Cov$returns
Pelet_Cov$convergence
c("Est_DF" = Pelet_Cov$Est_N, "Real_DF" = Tlength)
c("Test" = HypTestResHess$Test, "Sig Level" = HypTestResHess$Significance, "FWER Method" = HypTestResHess$MH_method)

HypTestResHess$Results[order(HypTestResHess$Results$Real),]
HypTestResGrad$Results[order(HypTestResGrad$Results$Real),]
HypTestResComb$Results[order(HypTestResComb$Results$Real),]
select(HypTestResComb$Results, Est., Lower, Upper) %>% linkFun$FUN()

wilksTest(Pelet_Cov, sampleData$healthy, sampleData$sick, linkFun = linkFun, dim_alpha = 1)

multiRes <- multipleComparison(healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                               Tlength = Pelet_Cov$Est_N, p.adjust.method = "BH") %>%
  round(3) %>% vector2triangle()

# corrplot(1 - multiRes, is.corr = F)
# corrplot(multiRes < 0.05, is.corr = F)
# select(HypTestResComb$Results, Est., 'Adj.P-val', Reject_H0)
tt.all <- Sys.time() - tt.all
# })