#Set parameters for Real Sigma
p <- 10
location.par <- 0
scale.par <- 1

#Set parameters for alpha
scale.alpha <- c(1,0.4)
plot(seq(0,1,0.01), dbeta(seq(0,1,0.01), scale.alpha[1], scale.alpha[2]), type = "l", xlab = "x",
     ylab = "CDF", main = paste("Beta CDF with parameters", scale.alpha[1],", ", scale.alpha[2]))

#Build Real Sigma
temp <- rcauchy(p, location.par, scale.par)
real.sigma <- temp%*%t(temp) + sqrt(mean(temp^2))*diag(p)
rm(temp, location.par, scale.par)

#Build Real Corr (Theta)
real.theta <- real.sigma / sqrt(diag(real.sigma)%*%t(diag(real.sigma)))

#Build Real Alpha
alpha <- rbeta(p, scale.alpha[1], scale.alpha[2])
create_alpha_mat <- function(VECT){
  pelet <- VECT%*%t(VECT)
  diag(pelet) <- rep(1, p)
  return(pelet)
}
alpha.mat <- create_alpha_mat(alpha)
rm(scale.alpha)

#Simulate Sample
create_correlation_matrices <- function(real_corr, noise, sample_size){
  pelet <- list()
  for(b in 1:sample_size){
    if(nrow(real_corr) != ncol(real_corr)) stop("Correlation matrix not p x p")
    m <- nrow(real_corr)
    
    cap <- min(noise*3.5, 0.5)
    noise.vect <- apply(cbind(rnorm((m*(m-1))/2,0,noise), rep(cap, (m*(m-1)))), 1, min)
    noise.vect <- apply(cbind(noise.vect, rep(-cap, (m*(m-1)))), 1, max)
    
    noise.mat <- matrix(ncol = m, nrow = m)
    count <- 1
    diag(noise.mat) <- rep(0, m)
    
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        noise.mat[i,j] <- noise.vect[count]
        noise.mat[j,i] <- noise.vect[count]
        count <- count + 1
      }
    }
    temp <- real_corr+noise.mat
    temp[which(temp>1)] <- 1
    pelet[[b]] <- real_corr+noise.mat
  }
  return(pelet)
}

healthy_N <- 4
sick_N <- 4

healthy <- create_correlation_matrices(real.theta, 0.01, healthy_N)
sick <- create_correlation_matrices(real.theta*alpha.mat, 0.01, sick_N)

#Calculate mean Correlation
calculate_mean_matrix <- function(matrix_list){
  temp <- matrix(0, ncol = ncol(matrix_list[[1]]), nrow = nrow(matrix_list[[1]]))
  returns <- length(matrix_list)
  for(i in 1:returns) temp <- temp + matrix_list[[i]]
  return(temp/returns)
}
triangle_to_vector <- function(MATR){
  if(nrow(MATR) != ncol(MATR)) stop("Matrix not p x p")
  return(as.vector(MATR[lower.tri(MATR)]))
}
vector_to_triangle <- function(VECT){
  m <- length(VECT)
  p <- 0.5*c(1+sqrt(1+8*m), 1-sqrt(1+8*m))
  p <- p[which( (p==round(p))&p==abs(p) )]
  if(length(p)==0) stop("Vect length does not fit size of triangular matrix")
  pelet <- matrix(0, ncol = p, nrow = p)
  pelet[lower.tri(pelet)] <- VECT
  pelet <- pelet+t(pelet)
  diag(pelet) <- 1
  return(pelet)
}

#Fisher Trasnformation
Fisher_Z <- function(r) 0.5*log((1+r)/(1-r))
plot(seq(-1,1,0.01), Fisher_Z(seq(-1,1,0.01)), type = "l")

sick_mean <- calculate_mean_matrix(sick)
healthy_mean <- calculate_mean_matrix(healthy)


differential.alpha <- function(Alpha,Theta) as.vector(t(Alpha)%*%(Theta*
                                                              (Theta*create_alpha_mat(Alpha)-sick_mean)*
                                                              (1-diag(p))))

for.optim <- function(ALPHA, THETA) sum(triangle_to_vector(THETA*create_alpha_mat(ALPHA)-sick_mean)^2)
trim.double <- function(X, lower, upper, add.noise = FALSE, noise.sd = 0.005){
  m <- length(X)
  pelet <- apply(rbind(apply(rbind(X, rep(upper,m)), 2, min), rep(lower,m)), 2, max)
  
  if(add.noise){
    trimmed <- which(pelet%%1 == 0)
    if(length(trimmed)>0) pelet[trimmed] <- pelet[trimmed]+rnorm(length(trimmed), 0, noise.sd)
  }
  return(pelet)
}


Steps <- list()
temp.theta <- healthy_mean
temp.alpha <- newton.raphson.multi(FUN = function(A) differential.alpha(A, temp.theta),
                                           X0 = rep(0.5,p),
                                           data = FALSE,
                                           persic = 6)
Steps[[1]] <- list(theta = temp.theta, alpha = trim.double(temp.alpha, 0.5, 1.2, TRUE))

limit <- c(MaxLoop = 500, Persic = 10^(-5))
i <- 1
distance <- 100
while((i<=limit["MaxLoop"])&(distance>limit["Persic"])){
  temp.theta <- (sick_N*sick_mean/create_alpha_mat(temp.alpha)+healthy_N*healthy_mean)/(sick_N+healthy_N)
  temp.alpha <- newton.raphson.multi(FUN = function(A) differential.alpha(A, temp.theta),
                                     X0 = trim.double(temp.alpha, 0.5, 1.5, TRUE),
                                     data = FALSE,
                                     persic = 6)
  Steps[[i+1]] <- list(theta = temp.theta, alpha = temp.alpha)
  distance <- sqrt(sum((Steps[[i+1]]$alpha-Steps[[i]]$alpha)^2))
  i <- i + 1
}
rbind(Steps[[i]]$alpha, alpha)
cbind(triangle_to_vector(real.theta),triangle_to_vector(Steps[[i]]$theta))
