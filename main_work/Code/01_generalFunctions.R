ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("abind", "corrplot", "data.table", "Matrix", "matrixcalc",
              "mvtnorm", "numDeriv", "parallel", "plotly", "profvis", "progress",
              "Rcpp", "R.matlab", "stats4", "tidyverse", "GGally")
ipak(packages)

promptForCores <- function(){
  newPrompt <- TRUE
  
  if(tolower(.Platform$OS.type) == "windows"){
    ncores <<- 1
    newPrompt <- FALSE
  }
  
  det <- detectCores()
  
  if("ncores" %in% objects(name = .GlobalEnv))
    if(ncores > 0 & ncores < det) newPrompt <- FALSE

  if(newPrompt){
    ncores <- 0
    userans <- "0"
  } else {
    userans <- "y"
  }
  while(!(userans %in% c("y", "Y"))){
    ncores <- 0
    while(ncores < 1 | ncores >= det){
      ncores <- readline(paste0(det, " cores where detected. Please enter number of cores to use: "))
      ncores <- floor(as.numeric(ncores))
    }
    if(ncores > 1){
      userans <- readline(paste0(det, " detected, ", ncores, " used. Confirm (y)? "))
    } else {
      userans <- "y"
    }
  }
  ncores <<- ncores
  message(paste0("R will use ", ncores, " cores. 'ncores' saved to global environemnt."))
}
ncores <- ifelse(tolower(.Platform$OS.type) == "windows", 1, detectCores() - 2)
promptForCores()

buildCL <- function(ncores = .GlobalEnv$ncores, packageList, dataList){
  J <- length(packageList)
  cl <<- makeCluster(ncores)
  
  for(j in 1:J){
    tmp <- packageList[j]
    clusterExport(cl = cl, "tmp", envir = environment())
    clusterEvalQ(cl=cl, library(tmp, character.only = T))
  }
  clusterExport(cl = cl, dataList, envir = environment())
  setDefaultCluster(cl = cl)
  message("Cluster opened, saved as 'cl' on global environment. Do not forget to stopCluster. \n")
  return(cl)
}

terminateCL <- function(silent = FALSE){
  if("cl" %in% objects(envir = .GlobalEnv)){
    stopCluster(cl)
    setDefaultCluster()
    rm(cl, envir =  .GlobalEnv)
    if(!silent) message("Cluster terminated.")
  } else {
    warning("Cluster not found and was not closed. If cluster is not saved as 'cl' you must close the cluster manually!")
  }
}

powerMatrix <- function(MATR, pow){
  eigenMat <- eigen(MATR)
  return(eigenMat$vectors %*% diag(eigenMat$values ^ pow) %*% t(eigenMat$vectors))
}

regularizeMatrix <- function(MATR, method = c("diag", "constant", "avg.diag", "increase.diag"), const = 1, OnlyIfSing = TRUE){
  method <- method[1]
  if(!is.square.matrix(MATR)) stop("Matrix is not square.")
  if(OnlyIfSing & !is.singular.matrix(MATR)) {
    message("Matrix is invertible and was not changed.")
    return(MATR)
  }
  
  p <- nrow(MATR)
  if(method == "diag") {
    pelet <- MATR + diag(p)
  }
  if(method == "constant"){
    if(const <= 0) stop("In method 'constant' const must be greater than 0.")
    pelet <- MATR + const*diag(p)
  }
  if(method == "avg.diag"){
    if(const < 0 | const > 1) stop("In method 'avg.diag' const must be between 0-1")
    pelet <- (1 - const)*MATR + const*mean(diag(MATR))*diag(p)
  } 
  if(method == "increase.diag"){
    if(const < 0 | const > 1) stop("In method 'avg.diag' const must be between 0-1")
    pelet <- (1 - const)*MATR + const*diag(diag(MATR))
  }
  
  if(is.singular.matrix(pelet)){
    warning("Matrix still not invertible.")
  } else{
    message("Matrix is singular and was regularized.")
  }
  return(pelet)
}

#Calculate mean Correlation
calculate_mean_matrix <- function(matrix_array, do.mean = TRUE){
  temp <- matrix(0, ncol = dim(matrix_array)[2], nrow = dim(matrix_array)[1])
  returns <- dim(matrix_array)[3]
  for(i in 1:returns) temp <- temp + matrix_array[,,i]
  if(do.mean) return(temp/returns)
  return(temp)
}

#Check stationarity/invertability of AR/MA process
checkInv <- function(coefs, perc = 0.001){
  polfun <- function(x) 1 - sum(coefs*x^(1:length(coefs)))
  x <- sign(sapply(seq(-1, 1, by = perc), polfun))
  return(all(x == 1) | all(x == -1))
}

#Generate weighted sum of matrices from array
summatrix <- function(ARRAY, index, constants, weights = FALSE){
  if(missing(constants)) constants <- rep(1, length(index))
  if(weights) constants <- constants/sum(constants)
  pelet <- matrix(0, nrow = dim(ARRAY)[1], ncol = dim(ARRAY)[2])
  for(i in 1:length(index)){
    pelet <- pelet + ARRAY[,,index[i]]*constants[i]
  }
  return(pelet)
}

#Generate weighted sum of vectors from matrices
sumvector <- function(MATR, index, constants, weights = FALSE){
  if(missing(constants)) constants <- rep(1, length(index))
  if(weights) constants <- constants/sum(constants)
  pelet <- numeric(ncol(MATR))
  for(i in 1:length(index)){
    pelet <- pelet + MATR[index[i],]*constants[i]
  }
  return(pelet)
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
cor.matrix_to_norm.matrix <- function(ARRAY) t(apply(ARRAY, 3, triangle2vector))

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
triangle2vector <- function(MATR , diag = FALSE){
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

vector2triangle <- function(VECT, diag = FALSE, truncdiag = 1){
  m <- length(VECT)
  
  one <- ifelse(diag, -1, 1)
  p <- 0.5*c(one + sqrt(1 + 8*m), one - sqrt(1 + 8*m))
  p <- p[which( (p==round(p)) & p==abs(p) )]
  if(length(p)==0) stop("Vect length does not fit size of triangular matrix")
  pelet <- matrix(0, ncol = p, nrow = p)
  pelet[lower.tri(pelet, diag = diag)] <- VECT
  
  if(diag){
    pelet <- pelet + t(pelet) - diag(diag(pelet))
  } else {
    pelet <- pelet + t(pelet)
    if(!is.null(truncdiag)) diag(pelet) <- truncdiag
  }
  return(pelet)
}

#Calculate Maholonobis norm of a vector. Default is regular norm.
vnorm <- function(x, MATR, sqroot = FALSE, solve_matr = FALSE){
  if(solve_matr) MATR <- solve(MATR)
  if(missing(MATR)) { pelet <- sum(x^2) } else{ pelet <- as.vector(t(x)%*%MATR%*%x) }
  if(sqroot) pelet <- sqrt(pelet)
  return(pelet)
}

#Create a symmetric matrix from a vector
# vector2triangle_old <- function(VECT){
#   m <- length(VECT)
#   p <- 0.5*c(1+sqrt(1+8*m), 1-sqrt(1+8*m))
#   p <- p[which( (p==round(p))&p==abs(p) )]
#   if(length(p)==0) stop("Vect length does not fit size of triangular matrix")
#   pelet <- matrix(0, ncol = p, nrow = p)
#   pelet[lower.tri(pelet)] <- VECT
#   pelet <- pelet+t(pelet)
#   diag(pelet) <- 1
#   return(pelet)
# }
# 
