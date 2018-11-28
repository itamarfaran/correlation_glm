source("Main Work/Code/generalFunctions.R")
source("Main Work/code/estimationFunctions2.R")
source("Main Work/code/simulationFunctions.R")

tt <- rep(Sys.time(), 2)
ncores <- detectCores() - 1 #### Enter here!
if(ncores > 1) requiredFunction <- c("Estimate.Loop", "Estimate.Loop2",
                                     "cor.matrix_to_norm.matrix", "triangle_to_vector","vector_to_triangle",
                                     "create_alpha_mat", "clean_sick", "vnorm", "compute_estimated_N","vector_var_matrix_calc_COR",
                                     "minusloglik", "bootstrapFunction")

ARMAdetails <- list(ARsick = c(0.4, -0.2), MAsick = c(0.4),
                    ARhealth = c(0.2, -0.1), MAhealth = c(0.4))
sapply(ARMAdetails, checkInv)

Tlength <- 115
B <- 100
p <- 10
sampleDataB <- createSamples(B = B, nH = 107, nS = 92, p = p, Tlength = Tlength,
                             percent_alpha = 0.4, range_alpha = c(0.6, 0.8),
                             ARsick = ARMAdetails$ARsick, MAsick = ARMAdetails$MAsick,
                             ARhealth = ARMAdetails$ARhealth, MAhealth = ARMAdetails$MAhealth)
