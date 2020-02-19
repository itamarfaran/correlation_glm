source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

sampleData <- prepare_corrmat_data(
  link = "main_work/Data/ADNI_data_AD_CN.mat",
  corr_matrix_name = 'all.corrmats',
  healthy_index_name = 'CONTROLS',
  sick_index_name = 'AD'
  )

test_corr_mat(sampleData)

linkFun <- linkFunctions$multiplicative_identity

tt.est <- Sys.time()
Pelet.Real <- estimateAlpha(healthy.data = sampleData$healthy, sick.data = sampleData$healthy,
                            T_thresh = 320, linkFun = linkFun, updateU = 1,
                            progress = TRUE)
tt.est <- Sys.time() - tt.est

save(Pelet.Real, file = paste0("main_work/Data/Enviroments/ADNIrun_", Sys.Date(), ".RData"))

tt.grad <- Sys.time()
fisherMatrGrad <- ComputeFisher(Pelet.Real, sampleData$sick, "Grad", linkFun = linkFun, ncores = ncores)  %>% regularizeMatrix()
tt.grad <- Sys.time() - tt.grad

save(Pelet.Real, fisherMatrGrad, file = paste0("main_work/Data/Enviroments/ADNIrun_", Sys.Date(), ".RData"))

tt.hess <- Sys.time()
fisherMatrHess <- ComputeFisher(Pelet.Real, sampleData$sick, "Hess", linkFun = linkFun)  %>% regularizeMatrix()
tt.hess <- Sys.time() - tt.hess

save(Pelet.Real, fisherMatrGrad, fisherMatrHess, file = paste0("main_work/Data/Enviroments/ADNIrun_", Sys.Date(), ".RData"))

fisherMatrComb <- fisherMatrHess %*% solve(fisherMatrGrad) %*% fisherMatrHess

HypTestResHess <- build_hyp.test(Pelet.Real, fisherMatrHess, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResGrad <- build_hyp.test(Pelet.Real, fisherMatrGrad, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResComb <- build_hyp.test(Pelet.Real, fisherMatrComb, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)

save(Pelet.Real, fisherMatrGrad, fisherMatrHess, fisherMatrComb,
     HypTestResHess, HypTestResGrad, HypTestResComb, file = paste0("main_work/Data/Enviroments/ADNIrun_", Sys.Date(), ".RData"))
