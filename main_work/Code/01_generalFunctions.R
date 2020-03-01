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


regularize_matrix <- function(MATR, method = c("diag", "constant", "avg.diag", "increase.diag"),
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
check_invertability_arma <- function(coefs, perc = 0.001){
  polfun <- function(x) 1 - sum(coefs*x^(1:length(coefs)))
  x <- sign(sapply(seq(-1, 1, by = perc), polfun))
  return(all(x == 1) | all(x == -1))
}
checkInv <- check_invertability_arma

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


mean_sqrt_diag <- function(x) mean(sqrt(diag(x)))
sqrt_diag <- function(x) sqrt(diag(x))


#Calculate non-biased estimates for Mean, Variance, Skewness and (Ex-)Kurtosis
central_moment <- function(x, norm = TRUE) {
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
corr_mat_array2normal_data_mat <- function(array) t(apply(array, 3, triangle2vector))


corr_mat_array2normal_data_mat_test <- function(obj, verbose = FALSE){
  if(class(obj) == "array"){
    message_ <- 'obj transformed from array to matrix'
    out <- corr_mat_array2normal_data_mat(obj)
  } else if (class(obj) %in% c("matrix", "data.frame")){
    message_ <- 'obj already in normal data matrix form'
    out <- obj
  } else {
    stop('obj not of class "array", "matrix" or "data.frame"')
  }
  if(verbose) message(message_)
  return(out)
}
  

cor.matrix_to_norm.matrix <- corr_mat_array2normal_data_mat


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


# Data preperation functions
prepare_corrmat_data <- function(link, corr_matrix_name, healthy_index_name, sick_index_name, subset, na_action = 'omit'){
  get_and_handle_na <- function(all_data, index, na_action){
    dta <- all_data[,,index]
    na_patients <- unique(which(is.na(dta), arr.ind = TRUE)[,3])
    if(length(na_patients)){
      if(na_action == 'omit'){
        dta <- dta[,,-na_patients]
      } else {
        stop('na_action must be one of \'omit\'')
      }
    }
    return(dta)
  }
  
  
  na_action <- na_action[1]
  real_dta <- readMat(link)
  corr_mats <- real_dta[[corr_matrix_name]]
  
  for(try in 1:10){
    if(class(corr_mats) == 'array') break()
    corr_mats <- simplify2array(corr_mats)
    if(try == 10) stop('corr_mats could not be simplified to array')
  }
  which_cols_na <- which(is.na(corr_mats[1,,1]), arr.ind = TRUE)
  p <- dim(corr_mats)[1] - length(which_cols_na)
  all_data <- array(dim = c(p, p, dim(corr_mats)[3]))
  if(length(which_cols_na)){
    for(i in 1:dim(all_data)[3]) all_data[,,i] <- force_symmetry(corr_mats[-which_cols_na, -which_cols_na, i])
  } else {
    for(i in 1:dim(all_data)[3]) all_data[,,i] <- force_symmetry(corr_mats[,,i])
  }
  
  healthy_dta <- get_and_handle_na(all_data, real_dta[[healthy_index_name]], na_action)
  sick_dta <- get_and_handle_na(all_data, real_dta[[sick_index_name]], na_action)
  
  if(!missing(subset)){
    healthy_dta <- healthy_dta[subset, subset, ]
    sick_dta <- sick_dta[subset, subset, ]
    p <- length(subset)
  }
  
  sample_data <- list(samples = list(healthy = healthy_dta, sick = sick_dta),
                      p = p, which_cols_na = which_cols_na, link = link)
  return(sample_data)
}


test_corr_mat <- function(dta){
  which_not_pos <- which(!apply(with(dta$samples, abind(healthy, sick)), 3, is.positive.definite))
  which_na <- list(
    healthy = which(is.na(dta$healthy), arr.ind = TRUE),
    sick = which(is.na(dta$sick), arr.ind = TRUE)
  )
  return(list(
    which_not_positive_definite = which_not_pos,
    which_na = which_na
  ))
}
