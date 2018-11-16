ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#Calculate mean Correlation
calculate_mean_matrix <- function(matrix_array){
  temp <- matrix(0, ncol = dim(matrix_array)[2], nrow = dim(matrix_array)[1])
  returns <- dim(matrix_array)[3]
  for(i in 1:returns) temp <- temp + matrix_array[,,i]
  return(temp/returns)
}

#Calculate non-biased estimates for Mean, Variance, Skewness and (Ex-)Kurtosis
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

#Take array of symmetric matrices and convert them to one data matrix
cor.matrix_to_norm.matrix <- function(ARRAY) t(apply(ARRAY, 3, triangle_to_vector))

#Build the alpha matrix according to the model
create_alpha_mat <- function(VECT){
  pelet <- VECT%*%t(VECT)
  diag(pelet) <- rep(1, length(VECT))
  return(pelet)
}

#Force Positive Definiteness
force_positive_definiteness <- function(MATR, sensitivity = 0.01, homoscedasticity = FALSE){
  if(!is.symmetric.matrix(MATR)) stop("MATR not symmetric")
  alpha_seq <- unique(c(0, seq(0,1, by = sensitivity), 1))
  pelet <- MATR
  if(homoscedasticity){
    if(mean(diag(MATR)) <= 0) stop("Diag mean not positive")
    diag_MATR <- mean(diag(MATR)) * diag(nrow(MATR))
  } else {
    if(any(diag(MATR) <= 0)) stop("Diag not positive")
    diag_MATR <- diag(diag(MATR))
  }
  
  i <- 1
  while(!is.positive.definite(pelet)){
    pelet <- alpha_seq[i]*diag_MATR + (1 - alpha_seq[i])*MATR
    i <- i+1
  }
  
  if(i == 1) return(list(Matrix = pelet,
                         Alpha = 0)) else return(list(Matrix = pelet,
                                                      Alpha = alpha_seq[i-1]))
}

#Force symmetry on non-symmetrical matrix
force_symmetry <- function(MATR) return((MATR + t(MATR))/2)

#Calculte sum of mahalonobis distances
SSS_norm.matrix <- function(DATA, mu, sigma, solve_sig = TRUE, reg.par = 0){
  if((reg.par < 0) | (reg.par > 1)) stop("reg.par not between [0,1]")
  sigma <- (1 - reg.par)*sigma + reg.par*mean(diag(sigma))*diag(length(mu))
  
  if(solve_sig) sigma <- solve(sigma)
  dist <- DATA - rep(1,nrow(DATA))%*%t(mu)
  return(sum(diag(dist%*%sigma%*%t(dist))))
}

#Retrieve lower/upper triangle of a matrix as a vector
triangle_to_vector <- function(MATR , diag = FALSE){
  if(nrow(MATR) != ncol(MATR)) stop("Matrix not p x p")
  return(as.vector(MATR[lower.tri(MATR, diag = diag)]))
}

#Trim extreme values
trim_num <- function(x, lower = -Inf, upper = Inf){
  pelet <- x
  pelet[x<lower] <- lower
  pelet[x>upper] <- upper
  return(pelet)
}

#Create a symmetric matrix from a vector
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

#Calculate Maholonobis norm of a vector. Default is regular norm.
vnorm <- function(x, MATR = NULL, sqroot = FALSE, solve_matr = FALSE){
  if(solve_matr) MATR <- solve(MATR)
  if(length(MATR)==0) { pelet <- sum(x^2) } else{ pelet <- as.vector(t(x)%*%MATR%*%x) }
  if(sqroot) pelet <- sqrt(pelet)
  return(pelet)
}

