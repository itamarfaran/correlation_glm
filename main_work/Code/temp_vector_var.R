# source("Main Work/Code/generalFunctions.R")
# source("Main Work/Code/estimationFunctions2.R")
# source("Main Work/Code/simulationFunctions.R")

p <- 50
MATR <- build_parameters(p, 0.5, c(0,1))$Corr.mat

vector_var_matrix_calc_COR2 <- function(MATR, cl = NULL, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  p <- dim(MATR)[1]
  m <- p*(p-1)/2
  tocomp <- triangle2vector(matrix(1:m^2, nrow = m, ncol = m), TRUE)
  order_vect <- cbind(unlist(lapply(1:(p - 1), function(i) rep(i, p - i))),
                      unlist(lapply(1:(p - 1), function(i) (i + 1):p)))
  
  MATR <- as.vector(MATR)
  
  # Need to work with triangle-vector here instead of matrix.. i.e moving from i,j to the correct
  # index in triangle2vector(MATR)
  real.cov2 <- function(i, j, k, l, MATR) {
    (MATR[(i-1)*p + j]*MATR[(k-1)*p + l]/2) * (MATR[(i-1)*p + k]^2 + MATR[(i-1)*p + l]^2 +
                                                 MATR[(j-1)*p + k]^2 + MATR[(j-1)*p + l]^2) -
      MATR[(i-1)*p + j]*(MATR[(i-1)*p + k]*MATR[(i-1)*p + l] + MATR[(j-1)*p + k]*MATR[(j-1)*p + l]) -
      MATR[(k-1)*p + l]*(MATR[(i-1)*p + k]*MATR[(j-1)*p + k] + MATR[(i-1)*p + l]*MATR[(j-1)*p + l]) +
      (MATR[(i-1)*p + k]*MATR[(j-1)*p + l] + MATR[(i-1)*p + l]*MATR[(j-1)*p + k])
  }
  
  fillmat <- function(k, MATR){
    indexes <- c(order_vect[ceiling(k/m),], order_vect[ifelse(k %% m == 0, m, k %% m),])
    return(real.cov2(indexes[1], indexes[2], indexes[3], indexes[4], MATR))
  }
  
  if(is.null(cl)){
    #pelet <- vector2triangle(unlist(lapply(tocomp, fillmat)), TRUE)
    pelet <- vector2triangle(sapply(tocomp, fillmat, MATR = MATR), TRUE)
    #pelet <- vector2triangle(vapply(tocomp, fillmat, FUN.VALUE = numeric(1) ), TRUE)
  } else {
    clusterExport(cl, c("m", "p", "order_vect", "real.cov2","MATR"), envir = environment())
    pelet <- vector2triangle(unlist(parLapply(cl, tocomp, fillmat)), TRUE)
  }
  
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

vector_var_matrix_calc_COR <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  real.cov2 <- function(i, j, k, l) {
    (MATR[i,j]*MATR[k,l]/2) * (MATR[i,k]^2 + MATR[i,l]^2 + MATR[j,k]^2 + MATR[j,l]^2) -
      MATR[i,j]*(MATR[i,k]*MATR[i,l] + MATR[j,k]*MATR[j,l]) -
      MATR[k,l]*(MATR[i,k]*MATR[j,k] + MATR[i,l]*MATR[j,l]) +
      (MATR[i,k]*MATR[j,l] + MATR[i,l]*MATR[j,k])
  }
  
  p <- dim(MATR)[1]
  m <- p*(p-1)/2
  
  order_vect <- cbind(unlist(lapply(1:(p - 1), function(i) rep(i, p - i))),
                      unlist(lapply(1:(p - 1), function(i) (i + 1):p)))
  
  pelet <- matrix(nrow = m, ncol = m)
  for(i in 1:m){
    for(j in i:m){
      indexes <- c(order_vect[i,], order_vect[j,])
      pelet[i,j] <- real.cov2(indexes[1], indexes[2], indexes[3], indexes[4])
      pelet[j,i] <- pelet[i,j]
    }
  }
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

tt1 <- Sys.time()
pelet <- vector_var_matrix_calc_COR(MATR)
tt1 <- Sys.time() - tt1

tt2 <- Sys.time()
pelet2 <- vector_var_matrix_calc_COR2(MATR)
tt2 <- Sys.time() - tt2

tt1
tt2

identical(round(pelet, 2), round(pelet2 ,2))
rm(pelet, pelet2)

# MATR[1,7]
# MATR[(1-1)*p + 7]
# 
# MATR[34,45]
# MATR[(34-1)*p + 45]
# MATR[i,45]
# MATR[(i-1)*p + j]

