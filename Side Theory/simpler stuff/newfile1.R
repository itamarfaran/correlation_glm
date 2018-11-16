library(MASS)
library(matrixcalc)

p <- 50
mu <- as.matrix(round(runif(p, -1, 1), 3))

n <- 100
X <- mvrnorm(n = n, mu = mu, Sigma = mu %*% t(mu) + diag(p))

minusloglike1 <- function(mu, DATA){
  n <- nrow(DATA)
  
  Sigma <- mu %*% t(mu) + diag(p)
  if(!is.positive.definite(Sigma)) stop("Sigma not positive definite")
  
  dist <- X - rep(1,n) %*% t(mu)
  log(det(Sigma)) + mean(diag(dist %*% solve(Sigma) %*% t(dist)))
}

minusloglike2 <- function(mu, DATA){
  mtr <- function(MATR) mean(diag(MATR))
  n <- nrow(DATA)
  dist <- X - rep(1,n)%*%t(mu)
  return(  log(1 + t(mu)%*%mu) +
             mtr(dist%*%t(dist)) -
             mtr(dist%*%mu%*%t(mu)%*%t(dist)) / (1 + t(mu)%*%mu) )
}

minusloglike1(mu, X)

minusloglike2(mu, X)

MLE_sol <- optim(par = rep(0, p), fn = function(M) minusloglike1(M, X), method = "L-BFGS-B")
plot(mu, MLE_sol$par)
abline(a = 0, b = 1)

bias <- MLE_sol$par - mu
plot(mu, bias)
abline(a = 0, b = 0)
summary(lm(bias ~ mu + I(mu^2)))

alpha <- seq(-0.99, 0.99, by = 0.01)
pelet <- matrix(nrow = length(alpha), ncol = 2)
pelet[,1] <- alpha
colnames(pelet) <- c("Value", "Likelihood")

index <- 1
for(i in 1:length(alpha)){
  pelet[i,2] <- -minusloglike1(replace(mu, index, alpha[i]), X)
}

plot(pelet[,1:2], type = "l", col = "black")
points(mu[index], -minusloglike1(mu, X), col = "blue")
points(mean(X[,index]), -minusloglike1(replace(mu, index, mean(X[,index])), X), col = "darkgreen")
points(MLE_sol$par[index],
       -minusloglike1(replace(mu, index, MLE_sol$par[index]), X), col = "red")

