source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

link <- "main_work/Data/NMDA_all_data_AAL90.mat"
Real.dta <- readMat(link)

#Which coloumns are NA? (Usually, 87-88)
which.drop <- which(is.na(Real.dta$group.all[1,,1]), arr.ind = TRUE)
p <- dim(Real.dta$group.all)[1] - length(which.drop)

All.data <- array(dim = c(p, p, dim(Real.dta$group.all)[3]))

#Clean data from unknown coloumns
for(i in 1:dim(All.data)[3]) All.data[,,i] <- force_symmetry(Real.dta$group.all[-which.drop,-which.drop,i])

#Some observations still have NAs. Remove those observations:
healthy.Real_t <- All.data[,,Real.dta$CONTROLS]
healthy.Real <- healthy.Real_t[,,-unique(which(is.na(healthy.Real_t), arr.ind = TRUE)[,3])]

sick.Real_t <- All.data[,,Real.dta$NMDA]
sick.Real <- sick.Real_t[,,-unique(which(is.na(sick.Real_t), arr.ind = TRUE)[,3])]

sampleData <- list(healthy = healthy.Real, sick = sick.Real, p = p, which.drop = which.drop, link = link)

rm(which.drop, p, All.data, healthy.Real_t, healthy.Real, sick.Real_t, sick.Real)

#Now, are all matrices good to go?
abind(sampleData$healthy, sampleData$sick) %>% apply(3, is.positive.definite) %>% all()
sum(is.na(sampleData$healthy)) + sum(is.na(sampleData$sick))

linkFun <- linkFunctions$multiplicative_identity

tt.est <- Sys.time()
Pelet.Real <- estimateAlpha(healthy.data = sampleData$healthy, sick.data = sampleData$healthy,
                            T_thresh = 320, linkFun = linkFun, updateU = 1,
                            progress = TRUE)
tt.est <- Sys.time() - tt.est

save(Pelet.Real, file = paste0("main_work/Data/Enviroments/NMDArun_", Sys.Date(), ".RData"))

tt.grad <- Sys.time()
fisherMatrGrad <- ComputeFisher(Pelet.Real, sampleData$sick, "Grad", linkFun = linkFun, ncores = ncores)  %>% regularizeMatrix()
tt.grad <- Sys.time() - tt.grad

save(Pelet.Real, fisherMatrGrad, file = paste0("main_work/Data/Enviroments/NMDArun_", Sys.Date(), ".RData"))

tt.hess <- Sys.time()
fisherMatrHess <- ComputeFisher(Pelet.Real, sampleData$sick, "Hess", linkFun = linkFun)  %>% regularizeMatrix()
tt.hess <- Sys.time() - tt.hess

save(Pelet.Real, fisherMatrGrad, fisherMatrHess, file = paste0("main_work/Data/Enviroments/NMDArun_", Sys.Date(), ".RData"))

fisherMatrComb <- fisherMatrHess %*% solve(fisherMatrGrad) %*% fisherMatrHess

HypTestResHess <- build_hyp.test(Pelet.Real, fisherMatrHess, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResGrad <- build_hyp.test(Pelet.Real, fisherMatrGrad, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)
HypTestResComb <- build_hyp.test(Pelet.Real, fisherMatrComb, linkFun = linkFun, sampleData$alpha, Real = sampleData$alpha)

save(Pelet.Real, fisherMatrGrad, fisherMatrHess, fisherMatrComb,
     HypTestResHess, HypTestResGrad, HypTestResComb, file = paste0("main_work/Data/Enviroments/NMDArun_", Sys.Date(), ".RData"))
