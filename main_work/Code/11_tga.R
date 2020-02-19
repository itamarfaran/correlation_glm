source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")


sampleData <- prepare_corrmat_data(
  link = "main_work/Data/Amnesia_all_AAL.mat",
  corr_matrix_name = 'corrmats',
  healthy_index_name = 'CONTROLS',
  sick_index_name = 'TGA'
)

test_corr_mat(sampleData)
linkFun <- linkFunctions$multiplicative_identity

Pelet.Real <- estimateAlpha(healthy.data = sampleData$healthy, sick.data = sampleData$healthy,
                            T_thresh = 320, linkFun = linkFun, updateU = 1,
                            progress = TRUE)

save(Pelet.Real, file = paste0("main_work/Data/Enviroments/TGArun_", Sys.Date(), ".RData"))

save(Pelet.Real, fisherMatrGrad, fisherMatrHess, fisherMatrComb,
     HypTestResHess, HypTestResGrad, HypTestResComb, file = paste0("main_work/Data/Enviroments/TGArun_", Sys.Date(), ".RData"))
