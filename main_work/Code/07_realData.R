source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions2.R")

link <- "Main Work/Data/NMDA_all_data_AAL90.mat"
Real.dta <- readMat(link)

#Which coloumns are NA? (Usually, 87-88)
which.drop <- which(is.na(Real.dta$group.all[1,,1]), arr.ind = TRUE)
p <- dim(Real.dta$group.all)[1] - length(which.drop)

All.data <- array(dim = c(p, p, dim(Real.dta$group.all)[3]))

#Clean data from unknown coloumns
for(i in 1:dim(All.data)[3]){
  All.data[,,i] <- force_symmetry(Real.dta$group.all[-which.drop,-which.drop,i])
}
dim(All.data)

#Some observations still have NAs. Remove those observations:
healthy.Real_t <- All.data[,,Real.dta$CONTROLS]
healthy.Real <- healthy.Real_t[,,-unique(which(is.na(healthy.Real_t), arr.ind = TRUE)[,3])]
rm(healthy.Real_t)

sick.Real_t <- All.data[,,Real.dta$NMDA]
sick.Real <- sick.Real_t[,,-unique(which(is.na(sick.Real_t), arr.ind = TRUE)[,3])]
rm(sick.Real_t)

#Now, are all matrices good to go?
abind(healthy.Real, sick.Real) %>% apply(3, is.positive.definite) %>% all()
sum(is.na(healthy.Real)) + sum(is.na(sick.Real))

Pelet.Real <- Estimate.Loop(healthy.Real, sick.Real)

Pelet.Real$theta
mean(Pelet.Real$alpha < 0.9)
