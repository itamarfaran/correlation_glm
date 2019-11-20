ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


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


packages <- c("abind", "corrplot", "data.table", "Matrix", "matrixcalc",
              "mvtnorm", "numDeriv", "parallel", "plotly", "profvis", "progress",
              "Rcpp", "R.matlab", "stats4", "tidyverse", "GGally", "pbmcapply", "pbapply")
ipak(packages)

ncores <- ifelse(tolower(.Platform$OS.type) == "windows", 1, detectCores() - 2)
promptForCores()


powerMatrix <- function(MATR, pow){
  eigenMat <- eigen(MATR)
  return(eigenMat$vectors %*% diag(eigenMat$values ^ pow) %*% t(eigenMat$vectors))
}


regularizeMatrix <- function(MATR, method = c("diag", "constant", "avg.diag", "increase.diag"),
                             const = 1, OnlyIfSing = TRUE){
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
calculate_mean_matrix <- function(matrix_array, do.mean = TRUE) summatrix(matrix_array, weights = do.mean)


#Check stationarity/invertability of AR/MA process
checkInv <- function(coefs, perc = 0.001){
  polfun <- function(x) 1 - sum(coefs*x^(1:length(coefs)))
  x <- sign(sapply(seq(-1, 1, by = perc), polfun))
  return(all(x == 1) | all(x == -1))
}


#Generate weighted sum of matrices from array
summatrix <- function(ARRAY, index = 1:(dim(ARRAY)[3]), constants = rep(1, length(index)), weights = FALSE){
  # if(missing(index)) index <- 1:(dim(ARRAY)[3])
  # if(missing(constants)) constants <- rep(1, length(index))
  if(weights) constants <- constants/sum(constants)
  pelet <- matrix(0, nrow = dim(ARRAY)[1], ncol = dim(ARRAY)[2])
  for(i in 1:length(index)){
    pelet <- pelet + ARRAY[,,index[i]]*constants[i]
  }
  return(pelet)
}


#Generate weighted sum of vectors from matrices
sumvector <- function(MATR, index = 1:(nrow(MATR)), constants = rep(1, length(index)), weights = FALSE){
  # if(missing(index)) index <- 1:(nrow(MATR))
  # if(missing(constants)) constants <- rep(1, length(index))
  if(weights) constants <- constants/sum(constants)
  pelet <- numeric(ncol(MATR))
  for(i in 1:length(index)){
    pelet <- pelet + MATR[index[i],]*constants[i]
  }
  return(pelet)
}


#Calculate non-biased estimates for Mean, Variance, Skewness and (Ex-)Kurtosis
central.moment <- function(x, norm = TRUE) {
  b <- numeric(4)
  names(b) <- c("Mean", "Variance", "Skewness",
                ifelse(norm, "Kurtosis", "Ex.Kurtosis"))

  n <- length(x)
  b[1] <- mean <- mean(x)
  b[2] <- var(x)
  sd <- sqrt(b[2] * ((n - 1)/n))
  norm_x <- (x - mean)/sd
  
  skew <- mean(norm_x^3)
  kurt <- mean(norm_x^4)
  
  b[3] <- (sqrt(n*(n - 1)) / (n - 2)) * skew
  b[4] <- ( (n - 1) / ((n - 2) * (n - 3)) ) * ((n + 1) * kurt + 6) + ifelse(norm, -3, 0)

  return(b)
}


#Take array of symmetric matrices and convert them to one data matrix
cor.matrix_to_norm.matrix <- function(ARRAY) t(apply(ARRAY, 3, triangle2vector))


#Build the alpha matrix according to the multiplicative model
create_alpha_mat <- function(alpha, dim_alpha = 1, diag_val = 1){
  alpha <- matrix(alpha, nc = dim_alpha)
  output <- alpha %*% t(alpha)
  diag(output) <- diag_val
  return(output)
}


create_additive_alpha_mat <- function(alpha, dim_alpha = 1, diag_val = 0){
  if(dim_alpha > 1) stop("currently support only dim_alpha = 1")
  output <- replicate(length(a), a) + t(replicate(length(a), a))
  diag(output) <- 0
  return(output)
}


#Force symmetry on non-symmetrical matrix
force_symmetry <- function(MATR) return((MATR + t(MATR))/2)


#Retrieve lower/upper triangle of a matrix as a vector
triangle2vector <- function(MATR , diag = FALSE){
  if(nrow(MATR) != ncol(MATR)) stop("Matrix not p x p")
  return(as.vector(MATR[lower.tri(MATR, diag = diag)]))
}


vector2triangle <- function(VECT, diag = FALSE, truncdiag = 1){
  m <- length(VECT)
  
  one <- ifelse(diag, -1, 1)
  p <- 0.5*c(one + sqrt(1 + 8*m), one - sqrt(1 + 8*m))
  p <- p[which( (p==round(p)) & p==abs(p) )]
  if(length(p)==0) stop("Vect length does not fit size of triangular matrix")
  
  output <- matrix(0, ncol = p, nrow = p)
  output[lower.tri(output, diag = diag)] <- VECT
  
  if(diag){
    output <- output + t(output) - diag(diag(output))
  } else {
    output <- output + t(output)
    if(!is.null(truncdiag)) diag(output) <- truncdiag
  }
  return(output)
}


#Calculate Maholonobis norm of a vector. Default is regular norm.
vnorm <- function(x, MATR, sqrt = FALSE, solve_matr = FALSE){
  if(solve_matr) MATR <- solve(MATR)
  if(missing(MATR)) { pelet <- sum(x^2) } else { pelet <- as.vector(t(x) %*% MATR %*% x) }
  if(sqrt) pelet <- sqrt(pelet)
  return(pelet)
}


#Calculte sum of mahalonobis distances
SSS_norm.matrix <- function(DATA, mu, sigma, solve_sig = TRUE, reg.par = 0){
  if((reg.par < 0) | (reg.par > 1)) stop("reg.par not between [0,1]")
  sigma <- (1 - reg.par)*sigma + reg.par*mean(diag(sigma))*diag(length(mu))
  
  if(solve_sig) sigma <- solve(sigma)
  dist <- DATA - rep(1,nrow(DATA))%*%t(mu)
  return(sum(diag(dist%*%sigma%*%t(dist))))
}

