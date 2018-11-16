central.moment <- function(x,norm=TRUE) {
  n<-length(x)
  b<-vector()
  
  mean<-mean(x)
  sd<-sqrt(var(x)*((n-1)/n))
  
  b<-c(b,mean)
  b<-c(b,var(x))
  
  skew<-mean((x-mean)^3)/(sd^3)
  kurt<-(mean((x-mean)^4)/(sd^4))
  
  b<-c(b,(sqrt(n*(n-1))/(n-2))*skew)
  b<-c(b,((n-1)/((n-2)*(n-3)))*((n+1)*kurt+6))
  if (norm) b[4]<-b[4]-3
  
  names(b) <- c("Mean", "Variance", "Skewness", "Kurtosis")
  if (norm) names(b)[4]<-"Ex.Kurtosis"
  return(b)
}

create_mles <- function(theta, n, B){
  mle_pelet <- numeric(length = B)
  suff_stat <- matrix(nrow = B, ncol = 4)
  colnames(suff_stat) <- c("meanX", "meanX^2", "sdX", "thetaHat")
  for(b in 1:B){
    tempsamp <- rnorm(n, theta, abs(theta))
    suff_stat[b,] <- c(mean(tempsamp), mean(tempsamp^2), sqrt(n/(n-1))*sd(tempsamp), 0) 
    tempmle <- -0.5*(suff_stat[b,1] + c(1, -1)*sqrt(suff_stat[b,1]^2 + 4*suff_stat[b,2]))
    suff_stat[b,4] <- tempmle[which.min((abs(tempmle) - suff_stat[b,3])^2 +
                                          (tempmle - suff_stat[b,1])^2)]
  }
  return(suff_stat)
}

theta <- runif(1, -8, 8)
n <- 22

theta_midgam <- create_mles(theta, n, 10000)
bias <- theta_midgam[,4]-theta
normal_bias <- (theta_midgam[,4]-theta)/sqrt((theta^2)/(3*n))

hist(normal_bias, col = "lightblue", probability = TRUE)
lines( density(normal_bias, adjust = 2), col ="darkgreen", lwd = 2)
central.moment(normal_bias)
kolmogorovtest(normal_bias)

colMeans(theta_midgam)
theta
abs(colMeans(abs(theta_midgam))[-2]-abs(theta))
which.min(abs(colMeans(abs(theta_midgam))[-2]-abs(theta)))

sd(theta_midgam[,4])
sd(theta_midgam[,1])
sd(theta_midgam[,3])
sqrt((theta^2)/(3*n))

sd(theta_midgam[,4])/sqrt((theta^2)/(3*n))