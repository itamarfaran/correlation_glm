source("Main Work/Code/generalFunctions.R")
source("Main Work/code/estimationFunctions2.R")
source("Main Work/code/simulationFunctions.R")

Tlength <- 115
sampleData <- createSamples(nH = 107, nS = 92, p = 12, Tlength = Tlength,
                                percent_alpha = 0.4, range_alpha = c(0.6, 0.8))

#Are all matrices positive definite?
all(abind(sampleData$healthy, sampleData$sick, along = 3) %>%
      apply(3, is.positive.definite))


Pelet_IID <- Estimate.Loop(sampleData$healthy, sampleData$sick, MaxLoop = 100)

# emp <- sampleData$healthy %>% cor.matrix_to_norm.matrix() %>% cov() %>% triangle_to_vector(diag = TRUE)
# theo <- sampleData$healthy %>% calculate_mean_matrix() %>% vector_var_matrix_calc_COR() %>%
#   triangle_to_vector(diag = TRUE)
# 
# data.frame(Theoretical = theo, Empirical = emp) %>% ggplot(aes(x = Theoretical, y = Empirical)) + 
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   geom_point(col = "blue", alpha = 0.5) + geom_smooth(method = "lm")

Pelet_Cov <- Estimate.Loop2(theta0 = Pelet_IID$theta, alpha0 = Pelet_IID$alpha,
                            healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                            T_thresh = Tlength, method = "Nelder-Mead")

fisherMatrHess <- ComputeFisher(Pelet_Cov, sampleData$sick, "Hess")
fisherMatrGrad <- ComputeFisher(Pelet_Cov, sampleData$sick, "Grad")

HypTestResHess <- build_hyp.test(Pelet_Cov, fisherMatrHess, sampleData$alpha, MH_method = "holm", const = 1, effectiveN = Pelet_Cov$Est_N)
HypTestResGrad <- build_hyp.test(Pelet_Cov, fisherMatrGrad, sampleData$alpha, MH_method = "holm", const = 1, effectiveN = Pelet_Cov$Est_N)

Pelet_Cov$returns
Pelet_Cov$convergence
c("Est_DF" = Pelet_Cov$Est_N, "Real_DF" = Tlength)
c("Test" = tmp$Test, "Sig Level" = tmp$Significance, "FWER Method" = tmp$MH_method)

HypTestResHess$Results[order(HypTestResHess$Results$Real),]
HypTestResGrad$Results[order(HypTestResGrad$Results$Real),]

wilksTest(Pelet_Cov, sampleData$healthy, sampleData$sick)



