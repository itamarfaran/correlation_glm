source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

link <- "main_work/Data/ADNI_data_AD_CN.mat"
Real.dta <- readMat(link)

corr_mats <- simplify2array(simplify2array(Real.dta$all.corrmats))
#Which coloumns are NA? (Usually, 87-88)
which.drop <- which(is.na(corr_mats[1,,1]), arr.ind = TRUE)
p <- dim(corr_mats)[1] - length(which.drop)

All.data <- array(dim = c(p, p, dim(corr_mats)[3]))
for(i in 1:dim(All.data)[3]) All.data[,,i] <- force_symmetry(corr_mats[,,i])

#Some observations still have NAs. Remove those observations:
healthy.Real_t <- All.data[,,Real.dta$CONTROLS]
healthy.Real <- healthy.Real_t[,,-unique(which(is.na(healthy.Real_t), arr.ind = TRUE)[,3])]

str(healthy.Real)

dim(healthy.Real)

cov_empirical <- healthy.Real %>% cor.matrix_to_norm.matrix() %>% cov() %>% triangle2vector(diag = F)
cov_theo <- healthy.Real %>% calculate_mean_matrix() %>% vector_var_matrix_calc_COR_C() %>% triangle2vector(diag = F)

sqrt(mean((cov_empirical - cov_theo)^2))
mean(abs(cov_empirical - cov_theo))

data.frame(emp = cov_empirical,
           theo = cov_theo) %>%
  sample_n(25000) %>%
  ggplot(aes(x = emp, y = theo)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1)

data.frame(
  emp = cov_empirical %>% vector2triangle(diag = T) %>% diag(),
  theo = cov_theo %>% vector2triangle(diag = T) %>% diag()
  ) %>%
  # sample_n(25000) %>%
  ggplot(aes(x = emp, y = theo)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1)




