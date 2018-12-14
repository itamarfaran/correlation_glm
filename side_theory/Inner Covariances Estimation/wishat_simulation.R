library(MASS)

##Build Sigma. Cauchy distribution is used to create extreme values.
location.par <- 0
scale.par <- 1
deg <- 8

x <- as.matrix(rcauchy(deg, location.par, scale.par))
real_sigma <- x%*%t(x)+diag(runif(deg, 0, max(x)))
rm(x, location.par, scale.par)

###Choose wishart's df and number of simulations.
wishart.df <- 10
returns <- 100

###Create list with matrix simulations.
results <- list()
for(n in 1:returns){
  x <- mvrnorm(1, rep(0, deg), real_sigma)
  temp.matrix <- x%*%t(x)
  for(j in 2:wishart.df){
    x <- mvrnorm(1, rep(0, deg), real_sigma)
    temp.matrix <- temp.matrix + x%*%t(x)
  }
  results[[n]] <- temp.matrix
}
rm(j, n, temp.matrix, x)


retrieve.object <- function(i, j, k, l, indexes = NULL){
  if(length(indexes)==0) indexes = 1:length(results)
  temp.function <- Vectorize(function(n, i, j) results[[n]][i,j])
  temp <- cbind(temp.function(indexes, i, j), temp.function(indexes, k, l))
  colnames(temp) <- c(paste(i,j), paste(k,l))
  return(temp)
}

real.cov <- function(i, j, k, l) wishart.df*(real_sigma[i,k]*real_sigma[j,l]
                                            + real_sigma[i,l]*real_sigma[j,k])

compare_estimate <- function(i, j, k, l){
  real <- real.cov(i, j, k, l)
  estimate <- cov(retrieve.object(i, j, k, l))[1,2]
  return(list(Index = c(i,j,k,l), Real = real, Estimate = estimate, Abs.dist = abs(real - estimate),
              Abs.dist.Real.Ratio = abs((real - estimate)/real)))
}

compare_estimate(1,1,1,1)

m <- 0.5*deg*(deg+1)
real.cov.matr <- matrix(nrow = 0.5*m*(m+1), ncol = 8)
colnames(real.cov.matr) <- c("i","j","k","l","Real.COV","Estimate.COV",
                             "Abs.Distance","Distance.Real.Ratio")
count <- 1
for(i in 1:deg){
  for(j in i:deg){
    k <- i
    for(l in j:deg){
      real.cov.matr[count,1:4] <- c(i,j,k,l)
      real.cov.matr[count,5] <- real.cov(i,j,k,l)
      real.cov.matr[count,6] <- cov(retrieve.object(i, j, k, l))[1,2]
      real.cov.matr[count,7] <- abs(real.cov.matr[count,5]-real.cov.matr[count,6])
      real.cov.matr[count,8] <- abs(real.cov.matr[count,7]/real.cov.matr[count,5])
      count <- count + 1
    }
    for(k in (i+1):deg){
      for(l in k:deg){
        real.cov.matr[count,1:4] <- c(i,j,k,l)
        real.cov.matr[count,5] <- real.cov(i,j,k,l)
        real.cov.matr[count,6] <- cov(retrieve.object(i, j, k, l))[1,2]
        real.cov.matr[count,7] <- abs(real.cov.matr[count,5]-real.cov.matr[count,6])
        real.cov.matr[count,8] <- abs(real.cov.matr[count,7]/real.cov.matr[count,5])
        count <- count + 1
        if(count%%10==0) cat(paste(",", count))
      }
    }
  }
}


compare_estimate(sample(1:deg, 1), sample(1:deg, 1), sample(1:deg, 1), sample(1:deg, 1))


mean(real.cov.matr[,8]<2)
x <- real.cov.matr[,8][real.cov.matr[,8]<2]
hist(x, breaks = 50, xlab = "(Real.par - Est.)/Real.par",
     main = "Normalized Absolute Distance of Covariance Parameters", probability = TRUE, col = "lightblue")
lines(density(x, adjust = 2.5), col = "red", lwd = 2)
abline(v = mean(x), col = "black", lwd = 2, lty = 2)
quantile(x, 0:10/10)
sqrt(mean(x^2)) #sqrt mean square error
mean(x) #mean absolute square error
