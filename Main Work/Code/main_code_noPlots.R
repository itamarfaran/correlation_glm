source("Main Work/Code/generalFunctions.R")
source("Main Work/code/estimationFunctions.R")
source("Main Work/code/simulationFunctions.R")

packages <- c("abind", "corrplot", "numDeriv", "matrixcalc", "R.matlab", "stats4", "tidyverse", "data.table")
ipak(packages)

#Set sample size
healthy_N <- 107
sick_N <- 92
Tlength <- 115
p <- 7 # max{p} = 24


parameters <- build_parameters(p, 0.4, c(0.85,0.95))
real.theta <- parameters$Corr.mat
real.sigma <- parameters$Cov.mat
alpha <- parameters$Alpha
alpha.mat <- create_alpha_mat(alpha)
rm(parameters)

#Simulate Sample

healthy <- create_correlation_matrices(real.theta, healthy_N, Tlength)
sick <- create_correlation_matrices(real.theta*alpha.mat, sick_N, Tlength)

#Are all matrices positive definite?
all(abind(healthy, sick, along = 3) %>%
      apply(3, is.positive.definite))


Pelet_IID <- Estimate.Loop(healthy, sick, MaxLoop = 100)

emp <- healthy %>% cor.matrix_to_norm.matrix() %>% cov() %>% triangle_to_vector(diag = TRUE)
theo <- healthy %>% calculate_mean_matrix() %>% vector_var_matrix_calc_COR() %>%
  triangle_to_vector(diag = TRUE)

temp_mod <- lm(theo ~ 0 + emp)
round(c("Coef" = temp_mod$coefficients, "Effective_N" = Tlength), 3)

emp2 <- healthy %>% cor.matrix_to_norm.matrix() %>% cov() %>% diag()
theo2 <- healthy %>% calculate_mean_matrix() %>% vector_var_matrix_calc_COR() %>% diag()

temp_mod_diag <- lm(theo2 ~ 0 + emp2)
round(c("Coef" = temp_mod_diag$coefficients, "Effective_N" = Tlength), 3)

data.frame(Theoretical = theo, Empirical = emp) %>% ggplot(aes(x = Theoretical, y = Empirical)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(col = "blue", alpha = 0.5)


Pelet_Cov <- Estimate.Loop2(theta0 = triangle_to_vector(Pelet_IID$Estimates$theta),
                            alpha0 = Pelet_IID$Estimates$alpha, Healthy.ARR = healthy,
                            Sick.ARR = sick, T_thresh = Tlength, method = "Nelder-Mead")


tmp <- build_hyp.test(Pelet_Cov, alpha, method = "holm", const = 2, effectiveN = Pelet_Cov$Est_N)

Pelet_Cov$returns
Pelet_Cov$convergence
c("Est_DF" = Pelet_Cov$Est_N, "Real_DF" = Tlength)

tmp$Results[order(tmp$Results$Real),]
tmp$DF
tmp$Test
tmp$Significance
tmp$method

healthy.data <- healthy %>% cor.matrix_to_norm.matrix()
sick.data <- sick %>% cor.matrix_to_norm.matrix()


#Do a wilks test (chi-square)
chisq <- -2*( minusloglik(theta = Pelet_Cov$theta,
                          alpha = Pelet_Cov$alpha,
                          healthy.data = healthy.data,
                          sick.data = sick.data) -
                minusloglik(theta = rbind(healthy.data, sick.data) %>% colMeans(),
                            alpha = rep(1, length(Pelet_Cov$alpha)),
                            healthy.data = healthy.data,
                            sick.data = sick.data) )

c("Chisq_val" = chisq, "DF" = length(Pelet_Cov$alpha),
  "Pval" = 1 - pchisq(chisq, length(Pelet_Cov$alpha)))

if(p <= 7){
  B <- 100
  simuldat <- list()
  alpha_simul <- matrix(nrow = B, ncol = p)
  
  for (b in 1:B){
    cat(paste0("\n \n b: ", b, ";"))
    healthy_tmp <- create_correlation_matrices(real.theta, healthy_N, Tlength)
    sick_tmp <- create_correlation_matrices(real.theta*alpha.mat, sick_N, Tlength)
    
    res_specified <- Estimate.Loop2(theta0 = triangle_to_vector(Pelet_IID$Estimates$theta),
                                    alpha0 = Pelet_IID$Estimates$alpha,
                                    Healthy.ARR = healthy_tmp, Sick.ARR = sick_tmp, T_thresh = Tlength,
                                    comp_hess = FALSE)
    
    simuldat[[b]] <- list(healthy = healthy_tmp, sick = sick_tmp, specified = res_specified)
    alpha_simul[b,] <- res_specified$alpha
    
  }
  
  Emp_vs_Theo <- data.frame(Theoretic = tmp$SD,
                            Empiric = sapply(1:p, function(i) sd(alpha_simul[,i]))) %>%
    mutate(Quotient = Empiric/Theoretic)
  
  find_covariate <- lm(Empiric ~ 0 + Theoretic, data = Emp_vs_Theo)
  
  hist(Emp_vs_Theo$Quotient, main = "Quotient Histogram", xlab = "Quotient", col = "lightblue")
  summary(Emp_vs_Theo$Quotient)
  
  
  ggplot(Emp_vs_Theo, aes(x = Theoretic, y = Empiric)) +
    geom_abline(slope = 1, intercept = 0, col = "darkblue", size = 1, linetype = 2) +
    geom_abline(slope = find_covariate$coefficients, intercept = 0,
                col = "darkred", size = 1, linetype = 2) + 
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
    geom_point(col = "black", size = 1.5) + 
    labs(title = "Empiric SD of alpha VS. Theoretical",
         x = "Theoretical SD", y = "Empiric SD")
  summary(find_covariate)
  
  Tlength
  Pelet_Cov$Est_N
  Pelet_Cov$Est_N/Tlength
}

link2 <- paste0("C:/Users/Itamar/Google Drive/Documents/#My Documents/Study/99 Other/Binyamini/Main Work/Data/Enviroments/enviroment ", "p", p, " ", "B", B, " ", gsub(":", "", Sys.time()), ".RData")
save.image(link2)
rm(link2)

p <- 5
parameters_new <- build_parameters(p, 0.4, c(0.85,0.95))
Tlist <- c(10, 30, 50, 70, 100, 120, 150, 170, 200, 250, 300, 400, 700, 1000, 1500, 2000, 3000, 4000)
lngth_Tlist <- length(Tlist)
B <- 100

alpha_simul <- array(dim = c(B, p, lngth_Tlist))
alpha_sd_est <- matrix(0, nrow = lngth_Tlist, ncol = p)
eff_n <- numeric(lngth_Tlist)

cnt <- 0
for(ti in Tlist){
  cnt <- cnt + 1
  cat("\n DF =", ti, "(", cnt, "/", lngth_Tlist, "):")
  
  for (b in 1:B){
    isb1 <- (b==1)
    healthy_tmp <- create_correlation_matrices(parameters_new$Corr.mat, healthy_N, ti)
    sick_tmp <- create_correlation_matrices(parameters_new$Corr.mat*create_alpha_mat(parameters_new$Alpha), sick_N, ti)
    
    res_specified <- Estimate.Loop2(theta0 = triangle_to_vector(parameters_new$Corr.mat),
                                    alpha0 = parameters_new$Alpha,
                                    Healthy.ARR = healthy_tmp, Sick.ARR = sick_tmp, T_thresh = ti*10,
                                    comp_hess = isb1, progress = FALSE)
    alpha_simul[b,,cnt] <- res_specified$alpha
    if(isb1){
      alpha_sd_est[cnt,] <- sqrt(diag(solve(res_specified$Hess)))
      eff_n[cnt] <- res_specified$Est_N
      
    }
    if((10*b/B)%%1 == 0) cat(100*b/B, "%, ")
  }
}

emp_sds <- matrix(nrow = lngth_Tlist, ncol = p)
coeffs <- numeric(lngth_Tlist)
for(i in 1:lngth_Tlist) {
  emp_sds[i, ] <- apply(alpha_simul[,,i], 2, sd)
  coeffs[i] <- lm(emp_sds[i,] ~ 0 + alpha_sd_est[i,])$coef
}

tmp <- numeric()
for(i in 1:lngth_Tlist) tmp <- c(tmp, rep(Tlist[i], p))
forplt <- data.frame(DF = factor(tmp), Empiric = as.vector(t(emp_sds)), Estimate = as.vector(t(alpha_sd_est)))
tmp <- rbind(eff_n, round(eff_n/Tlist, 3), coeffs)
colnames(tmp) <- Tlist
row.names(tmp) <- c("Est_n", "Ratio", "Coeffs")

Show the Results
```{r}
ggplot(forplt, aes(x = Estimate, y = Empiric, col = DF)) + geom_point() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, size = 0.8, linetype = 2, col = "darkgrey") + 
  #geom_smooth(method = "lm", formula = y ~ 0 + x ,se = FALSE, linetype = 4) + 
  geom_abline(slope = 2, intercept = 0, size = 0.8, col = "darkred", linetype = 2)


data.frame(DF = Tlist, Coefficient = coeffs) %>% ggplot(aes(x = DF, y = Coefficient)) + geom_point() +
  geom_hline(yintercept = 1, col = "darkgrey") + geom_vline(xintercept = 0) + 
  geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE)

#View(t(tmp))
x <- 1/Tlist
summary(lm(coeffs ~ x))

df_error <- tmp[2,]-1
hist(df_error)
summary(df_error)
plot(Tlist, df_error)
summary(lm(df_error ~ 0 + Tlist))

link2 <- paste0("C:/Users/Itamar/Google Drive/Documents/#My Documents/Study/99 Other/Binyamini/Main Work/Data/Enviroments/enviroment robustness_check p_",
                p, " T_", lngth_Tlist, " ", gsub(":", "", Sys.time()), ".RData")
save.image(link2)
rm(link2)

