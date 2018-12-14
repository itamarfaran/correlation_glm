library(abind)
library(dplyr)
library(ggplot2)
library(matrixcalc)
library(tidyr)
library(MASS)

force_symmetry <- function(MATR) return((MATR + t(MATR))/2)

##Build Sigma. Cauchy distribution is used to create extreme values.
location.par <- 0
scale.par <- 1
deg <- 8

x <- as.matrix(rcauchy(deg, location.par, scale.par))
real_sigma <- x%*%t(x)+diag(runif(deg, 0, max(x)))
rm(x, location.par, scale.par)

is.positive.definite(real_sigma)

real_corr <- force_symmetry(cov2cor(real_sigma))

is.positive.definite(real_corr)

###Choose wishart's df and number of simulations.
wishart.df <- 10
returns <- 1000

###Create list with matrix simulations.
results <- list()
for(n in 1:returns){
  x <- mvrnorm(1, rep(0, deg), real_sigma)
  temp.matrix <- x%*%t(x)
  for(j in 2:wishart.df){
    x <- mvrnorm(1, rep(0, deg), real_sigma)
    temp.matrix <- temp.matrix + x%*%t(x)
  }
  results[[n]] <- temp.matrix/sqrt(diag(temp.matrix)%*%t(diag(temp.matrix)))
}
rm(j, n, temp.matrix, x)

retrieve.object <- function(i, j, k, l, indexes = NULL){
  if(length(indexes)==0) indexes = 1:length(results)
  temp.function <- Vectorize(function(n, i, j) results[[n]][i,j])
  return(cbind(temp.function(indexes, i, j), temp.function(indexes, k, l)))
}

real.cov <- function(i, j, k, l) (real_corr[i,k]*real_corr[j,l]
                                             + real_corr[i,l]*real_corr[j,k])/wishart.df

m <- 0.5*deg*(deg+1)
real.cov.matr <- matrix(nrow = 0.5*m*(m+1), ncol = 8)
colnames(real.cov.matr) <- c("i","j","k","l","Real.COV","Estimate.COV", "Distance","Distance.Real.Ratio")
count <- 1
for(i in 1:(deg-1)){
  for(j in (i+1):deg){
    k <- i
    for(l in j:deg){
      real.cov.matr[count,1:4] <- c(i,j,k,l)
      real.cov.matr[count,5] <- real.cov(i,j,k,l)
      real.cov.matr[count,6] <- cov(retrieve.object(i, j, k, l))[1,2]
      real.cov.matr[count,7] <- abs(real.cov.matr[count,5]-real.cov.matr[count,6])
      real.cov.matr[count,8] <- abs(real.cov.matr[count,7]/real.cov.matr[count,5])
      count <- count + 1
    }
    for(k in (i+1):(deg-1)){
      for(l in (k+1):deg){
        real.cov.matr[count,1:4] <- c(i,j,k,l)
        real.cov.matr[count,5] <- real.cov(i,j,k,l)
        real.cov.matr[count,6] <- cov(retrieve.object(i, j, k, l))[1,2]
        real.cov.matr[count,7] <- real.cov.matr[count,5]-real.cov.matr[count,6]
        real.cov.matr[count,8] <- real.cov.matr[count,7]/abs(real.cov.matr[count,5])
        count <- count + 1
        if(count%%10==0) cat(paste(",", count))
      }
    }
  }
}

real.cov.matr <- real.cov.matr[1:sum(!is.na(real.cov.matr[,5])), ]

compare_estimate <- function(i, j, k, l){
  real <- real.cov(i, j, k, l)
  estimate <- cov(retrieve.object(i, j, k, l))[1,2]
  return(list(Index = c(i,j,k,l), Real = real, Estimate = estimate, Abs.dist = abs(real - estimate),
              Abs.dist.Real.Ratio = abs((real - estimate)/real)))
}



compare_estimate(sample(1:deg, 1), sample(1:deg, 1), sample(1:deg, 1), sample(1:deg, 1))
compare_estimate(1,2,1,2)
compare_estimate(1,2,1,3)

hist(real.cov.matr[,8], breaks = 50, xlab = "(Real.par - Est.)/Real.par",
     main = "Normalized Absolute Distance of Correlation Parameters", probability = TRUE, col = "purple")
lines(density(real.cov.matr[,8], adjust = 2.5), col = "red", lwd = 2)
abline(v = mean(real.cov.matr[,8]), col = "black", lwd = 2, lty = 2)
quantile(real.cov.matr[,8], 0:10/10)
sqrt(mean(real.cov.matr[,7]^2)) #sqrt mean square error
mean(real.cov.matr[,7]) #mean absolute square error

variance.matrix <- matrix(0, ncol = deg, nrow = deg)
for(i in 1:(deg-1)){
  for(j in (i+1):deg) variance.matrix[i,j] <- cov(retrieve.object(i, j, i, j))[1,2]
}

variance.matrix <- variance.matrix + t(variance.matrix)
diag(variance.matrix) <- 1

distance.from.paramater <- sapply(results, function(X) triangle_to_vector((X-real_corr)/variance.matrix))
distance.from.paramater
hist(distance.from.paramater[,1])
