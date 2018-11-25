source("Main Work/Code/generalFunctions.R")
source("Main Work/code/estimationFunctions.R")
source("Main Work/code/simulationFunctions.R")
source("Main Work/code/games_with_QMLE.R")

#Set sample size
healthy_N <- 107
sick_N <- 92
Tlength <- 115
p <- 7


parameters <- build_parameters(p, 0.4, c(0.85,0.95))
real.theta <- parameters$Corr.mat
real.sigma <- parameters$Cov.mat
alpha <- parameters$Alpha
alpha.mat <- create_alpha_mat(alpha)
rm(parameters)

B <- 100
simuldat <- list()
alpha_simul <- matrix(nrow = B, ncol = p)

for (b in 1:B){
  cat(paste0("\n \n b: ", b, ";"))
  healthy_tmp <- create_correlation_matrices(real.theta, healthy_N, Tlength)
  sick_tmp <- create_correlation_matrices(real.theta*alpha.mat, sick_N, Tlength)
  
  res_unspecified <- Estimate.Loop(Healthy_List = healthy_tmp, Sick_List = sick_tmp)
  
  res_specified <- Estimate.Loop2(theta0 = triangle_to_vector(res_unspecified$Estimates$theta),
                                  alpha0 = res_unspecified$Estimates$alpha,
                                  Healthy.ARR = healthy_tmp, Sick.ARR = sick_tmp, T_thresh = Tlength,
                                  comp_hess = (b == 1))
  
  if(b == 1){
    b1Hess <- res_specified$Hess
    Bmatrix <- computeBmatr(sick_tmp, res_specified)
  }
  simuldat[[b]] <- list(healthy = healthy_tmp, sick = sick_tmp, specified = res_specified)
  alpha_simul[b,] <- res_specified$alpha
  
}

Emp_vs_Theo <- data.frame(Theoretic = sqrt(diag(solve(b1Hess))),
                          Empiric = sapply(1:p, function(i) sd(alpha_simul[,i]))) %>%
  mutate(Quotient = Empiric/Theoretic)

Emp_vs_Theo2 <- data.frame(Theoretic = sqrt(diag(solve(Bmatrix))),
                          Empiric = sapply(1:p, function(i) sd(alpha_simul[,i]))) %>%
  mutate(Quotient = Empiric/Theoretic)

grad_vs_hess <- data.frame(Grad = triangle_to_vector(solve(Bmatrix)),
                           Hess = triangle_to_vector(solve(b1Hess))) %>% 
  mutate(Quotient = Grad/Hess)

plot(grad_vs_hess[,1:2])
grad_vs_hess$Quotient
summary(lm(Grad ~ 0 + Hess, data = grad_vs_hess))



find_covariate <- lm(Empiric ~ 0 + Theoretic, data = Emp_vs_Theo2)

hist(Emp_vs_Theo2$Quotient, main = "Quotient Histogram", xlab = "Quotient", col = "lightblue")
summary(Emp_vs_Theo2$Quotient)

ggplot(Emp_vs_Theo2, aes(x = Theoretic, y = Empiric)) +
  geom_abline(slope = 1, intercept = 0, col = "darkblue", size = 1, linetype = 2) +
  geom_abline(slope = find_covariate$coefficients, intercept = 0,
              col = "darkred", size = 1, linetype = 2) + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  geom_point(col = "black", size = 1.5) + 
  labs(title = "Empiric SD of alpha VS. Theoretical",
       x = "Theoretical SD", y = "Empiric SD")
summary(find_covariate)

Tlength
res_specified$Est_N
res_specified$Est_N/Tlength

link2 <- paste0("Main Work/Data/Enviroments/enviroment ",
                "p", p, " ", "B", B, " ", gsub(":", "", Sys.time()), ".RData")
save.image(link2)
rm(link2)

p <- 5
parameters_new <- build_parameters(p, 0.4, c(0.85,0.95))
Tlist <- c(10, 30, 50, 70, 100, 120, 150, 170, 200, 250, 300, 400, 700, 1000, 1500, 2000, 3000, 4000)
lngth_Tlist <- length(Tlist)
B <- 100

alpha_simul <- array(dim = c(B, p, lngth_Tlist))
alpha_sd_est <- matrix(0, nrow = lngth_Tlist, ncol = p)
alpha_sd_est2 <- matrix(0, nrow = lngth_Tlist, ncol = p)
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
      alpha_sd_est2[cnt,] <- sqrt(diag(solve(computeBmatr(sick_tmp, res_specified))))
      eff_n[cnt] <- res_specified$Est_N
      
    }
    if((10*b/B)%%1 == 0) cat(100*b/B, "%, ")
  }
}

emp_sds <- matrix(nrow = lngth_Tlist, ncol = p)
coeffs <- numeric(lngth_Tlist)
for(i in 1:lngth_Tlist) {
  emp_sds[i, ] <- apply(alpha_simul[,,i], 2, sd)
  coeffs[i] <- lm(emp_sds[i,] ~ 0 + alpha_sd_est2[i,])$coef
}

sdsDims <- dim(emp_sds)
forLm <- matrix(0, nrow = prod(sdsDims), ncol = 4)
for(i in 1:sdsDims[1]){
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 1] <- Tlist[i]
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 2] <- alpha_sd_est[i,]
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 3] <- alpha_sd_est2[i,]
  forLm[((i-1)*sdsDims[2] + 1): (i*sdsDims[2]), 4] <- emp_sds[i,]
}

forLm <- as.data.frame(forLm)
colnames(forLm) <- c("DF", "EstimatedHess", "EstimatedGrad", "Empiric")
forLm <- mutate(forLm, QuotentHess = EstimatedHess/Empiric,
                QuotentGrad = EstimatedGrad/Empiric,
                QuotentBetween = EstimatedGrad/EstimatedHess)

summary(aov(QuotentGrad ~ factor(DF), data = forLm))
summary(aov(QuotentBetween ~ factor(DF), data = forLm))
summary(lm(QuotentBetween ~ log(DF), data = forLm))

summary(lm(QuotentGrad ~ 0 + log(DF)*QuotentHess, data = forLm))


ggplot(forLm, aes(x = factor(DF), col = factor(DF), y = QuotentGrad)) + geom_point()
ggplot(forLm, aes(x = factor(DF), col = factor(DF), y = QuotentHess)) + geom_point()
ggplot(forLm, aes(x = factor(DF), col = factor(DF), y = QuotentBetween)) + geom_point()

tmp <- numeric()
for(i in 1:lngth_Tlist) tmp <- c(tmp, rep(Tlist[i], p))
forplt <- data.frame(DF = factor(tmp), Empiric = as.vector(t(emp_sds)), Estimate = as.vector(t(alpha_sd_est2)))
tmp <- rbind(eff_n, round(eff_n/Tlist, 3), coeffs)
colnames(tmp) <- Tlist
row.names(tmp) <- c("Est_n", "Ratio", "Coeffs")


ggplot(forplt, aes(x = Estimate, y = Empiric)) + geom_point(aes(col = DF)) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(slope = 1, intercept = 0, size = 0.8, linetype = 2, col = "darkgrey") +
  geom_smooth(method = "lm", formula = y ~ 0 + x + exp(x) ,se = FALSE, linetype = 4)
  # geom_abline(slope = 2, intercept = 0, size = 0.8, col = "darkred", linetype = 2)


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

link2 <- paste0("Main Work/Data/Enviroments/enviroment robustness_check p_",
                p, " T_", lngth_Tlist, " ", gsub(":", "", Sys.time()), ".RData")
save.image(link2)
rm(link2)

