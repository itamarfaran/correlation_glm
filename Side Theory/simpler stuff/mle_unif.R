library(ggplot2)
library(tidyr)
library(dplyr)

f1 <- function(n, theta, epsilon) (1 - epsilon/theta)^n
f2 <- function(n, theta, epsilon) 2*pnorm((-1)*(epsilon * sqrt(3*n))/theta)

f1(4, 1, 0.5)
f2(4, 1, 0.5)

create_lines <- function(N, theta, epsilon){
  n <- length(N)
  pelet <- matrix(nrow = n, ncol = 3)
  colnames(pelet) <- c("n", "f1", "f2")
  
  pelet[,1] <- N
  pelet[,2] <- f1(n = N, theta = theta, epsilon = epsilon)
  pelet[,3] <- f2(n = N, theta = theta, epsilon = epsilon)
  
  return(pelet)
}

draw_lines <- function(MATR){
  MATR2 <- as.data.frame(MATR)
  colnames(MATR2) <- c("n", "Max{xi}", "Mean{X}")
  gather(MATR2, key = "Estimator", value = "Prob", -n) %>%
    ggplot(aes(x = n, y = Prob, color = Estimator)) + geom_line() + 
    geom_hline(yintercept = 0, size = 1) + geom_hline(yintercept = 1, color = "grey") + 
    labs(title = "Convergence of Estimators", y = "P(|Est. - Par.| > eps)")
}

seq(1, 50, by = 1) %>% create_lines(theta = 1, epsilon = 0.1) %>% draw_lines

seq(1, 300, by = 1) %>% create_lines(theta = 1, epsilon = 0.1) %>% draw_lines
