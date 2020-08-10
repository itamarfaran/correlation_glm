ipak <- function(..., only_install = FALSE){
  pkg <- unlist(list(...))
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  if(!only_install) sapply(pkg, require, character.only = TRUE)
}

packages <- c("abind", "corrplot", "data.table", "Matrix", "matrixcalc",
              "mvtnorm", "numDeriv", "parallel", "plotly", "profvis", "progress",
              "Rcpp", "R.matlab", "stats4", "tidyverse", "GGally", "pbmcapply", "pbapply")
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


matrix_pow <- function(x, pow){
  out <- with(eigen(x), vectors %*% diag(values ^ pow) %*% t(vectors))
  return(out)
}


regularize_matrix <- function(matr, method = c("constant", "avg_diag", "increase_diag"),
                             const = 1, only_if_singular = TRUE){
  method <- match.arg(method, c("constant", "avg_diag", "increase_diag"))
  if(!is.square.matrix(matr)) stop("Matrix is not square.")
  if(only_if_singular & !is.singular.matrix(matr)) {
    message("Matrix is invertible and was not changed.")
    return(matr)
  }
  
  p <- nrow(matr)
  if(method == "constant"){
    if(const <= 0) stop("In method 'constant' const must be greater than 0.")
    out <- matr + const*diag(p)
  }
  if(method == "avg_diag"){
    if(const < 0 | const > 1) stop("In method 'avg_diag' const must be between 0-1")
    out <- (1 - const)*matr + const*mean(diag(matr))*diag(p)
  } 
  if(method == "increase_diag"){
    if(const < 0 | const > 1) stop("In method 'increase_diag' const must be between 0-1")
    out <- (1 - const)*matr + const*diag(diag(matr))
  }
  
  if(is.singular.matrix(out)){
    warning("Matrix still not invertible.")
  } else{
    if(only_if_singular) message("Matrix is singular and was regularized.")
  }
  return(out)
}



#Check stationarity/invertability of AR/MA process
check_invertability_arma <- function(coefs, perc = 1e-03){
  polfun <- function(x) 1 - sum(coefs*x^(seq_along(coefs)))
  x <- sign(sapply(seq(-1, 1, by = perc), polfun))
  return(all(x == 1) | all(x == -1))
}


#Generate weighted sum of matrices from array
matrix_sum <- function(array_, index = 1:(dim(array_)[3]), constants = rep(1, length(index)), weights = FALSE){
  out <- matrix(0, nrow = dim(array_)[1], ncol = dim(array_)[2])
  
  if(weights) constants <- constants/sum(constants)
  
  for(i in seq_along(index)) out <- out + constants[i]*array_[,, index[i]]
  
  return(out)
}

#Calculate mean Correlation
calculate_mean_matrix <- function(matrix_array) matrix_sum(matrix_array, weights = TRUE)

#Generate weighted sum of vectors from matrices
vector_sum <- function(matr, index = 1:(nrow(matr)), constants = rep(1, length(index)), weights = FALSE, by_row = TRUE){
  out <- numeric(ncol(matr))
  
  if(weights) constants <- constants/sum(constants)
  
  for(i in seq_along(index)) out <- out + constants[i]*matr[index[i],]
  return(out)
}


mean_sqrt_diag <- function(x) mean(sqrt(diag(x)))
sqrt_diag <- function(x) sqrt(diag(x))


#Calculate non-biased estimates for Mean, Variance, Skewness and (Ex-)Kurtosis
central_moment <- function(x, norm = TRUE) {
  out <- numeric(4)
  names(out) <- c("Mean", "Variance", "Skewness", ifelse(norm, "Kurtosis", "Ex_Kurtosis"))

  n <- length(x)
  out[1] <- m <- mean(x)
  out[2] <- v <- var(x)
  
  s <- sqrt(v * (1 - 1/n))
  norm_x <- (x - m)/s
  
  skew <- mean(norm_x^3)
  kurt <- mean(norm_x^4)
  
  out[3] <- (sqrt(n*(n - 1)) / (n - 2)) * skew
  out[4] <- ( (n - 1) / ((n - 2) * (n - 3)) ) * ((n + 1) * kurt + 6) - 3*norm

  return(out)
}


#Take array of symmetric matrices and convert them to one data matrix
convert_corr_array_to_data_matrix <- function(array_) t(apply(array_, 3, triangle2vector))


convert_corr_array_to_data_matrix_test <- function(obj, verbose = FALSE){
  if(class(obj) == "array"){
    message_ <- 'obj transformed from array to matrix'
    out <- convert_corr_array_to_data_matrix(obj)
  } else if (class(obj) %in% c("matrix", "data.frame")){
    message_ <- 'obj already in normal data matrix form'
    out <- obj
  } else {
    stop('obj not of class "array", "matrix" or "data.frame"')
  }
  if(verbose) message(message_)
  return(out)
}
  

#Force symmetry on non-symmetrical matrix
force_symmetry <- function(matr) return((matr + t(matr))/2)


#Retrieve lower/upper triangle of a matrix as a vector
triangle2vector <- function(matr, diag = FALSE){
  if(nrow(matr) != ncol(matr)) stop("Matrix not p x p")
  return(as.vector(matr[lower.tri(matr, diag = diag)]))
}


vector2triangle <- function(vect, diag = FALSE, diag_value = NA){
  m <- length(vect)
  
  one <- ifelse(diag, -1, 1)
  p <- 0.5*c(one + sqrt(1 + 8*m), one - sqrt(1 + 8*m))
  p <- p[which( (p == round(p)) & p == abs(p) )]
  if(length(p) == 0) stop("Vect length does not fit size of triangular matrix")
  
  out <- matrix(0, ncol = p, nrow = p)
  out[lower.tri(out, diag = diag)] <- vect
  
  out <- out + t(out) - diag(diag(out))
  if(!diag) diag(out) <- diag_value
  return(out)
}


#Calculate Maholonobis norm of a vector. Default is regular norm.
vnorm <- function(x, matr, sqrt = FALSE, solve_matr = FALSE){
  if(missing(matr)){
    out <- sum(x^2)
  } else {
    if(solve_matr) matr <- solve(matr)
    out <- as.vector(t(x) %*% matr %*% x)
  }
  if(sqrt) out <- sqrt(out)
  return(out)
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
  real_dta <- R.matlab::readMat(link)
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

cov_known_mu <- function(x, y = NULL, mu, mu_x, mu_y, na.rm = FALSE, use = "everything"){

  if(is.vector(x)){
    if(na.rm) x <- x[!is.na(x)]
    if(is.null(y)) {
      y <- x
      if(missing(mu_x)) mu_x <- mu
      mu_y <- mu_x
    }
    out <- mean((x - mu_x) * (y - mu_y))
  }
  else if (is.matrix(x) | is.data.frame(x)){
    stop('var_known_mu supports univariate') # currently univariate
  } else {
    stop('x must be a vector, matrix a data.frame')
  }
  
  return(out)
}

var_known_mu <- function(x, mu = 0, na.rm = FALSE, use = "everything")
  cov_known_mu(x = x, y = NULL, mu_x = mu, na.rm = na.rm, use = use)

sd_known_mu <- function(x, mu, na.rm = FALSE)
  sqrt(var_known_mu(if (is.vector(x) || is.factor(x)) x else as.double(x), mu = mu, na.rm = na.rm))
