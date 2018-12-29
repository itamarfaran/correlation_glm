source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions2.R")

p <- 50
MATR <- build_parameters(p, 0.5, c(0,1))$Corr.mat

vector_var_matrix_calc_COR <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  real.cov2 <- function(i, j, k, l, MATR) {
    MATRij <- MATR[i,j]
    MATRkl <- MATR[k,l]
    MATRik <- MATR[i,k]
    MATRil <- MATR[i,l]
    MATRjk <- MATR[j,k]
    MATRjl <- MATR[j,l]
    
    (MATRij*MATRkl/2) * (MATRik^2 + MATRil^2 + MATRjk^2 + MATRjl^2) -
      MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
      MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
      (MATRik*MATRjl + MATRil*MATRjk)
  }
  
  p <- dim(MATR)[1]
  m <- p*(p-1)/2
  
  order_vect <- cbind(unlist(lapply(1:(p - 1), function(i) rep(i, p - i))),
                      unlist(lapply(1:(p - 1), function(i) (i + 1):p)))
  
  pelet <- matrix(nrow = m, ncol = m)
  for(i in 1:m){
    for(j in i:m){
      indexes <- c(order_vect[i,], order_vect[j,])
      pelet[i,j] <- real.cov2(indexes[1], indexes[2], indexes[3], indexes[4], MATR)
      pelet[j,i] <- pelet[i,j] # todo : check for profiler which tells you how much takes every line
    }
  }
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
} # todo : try exporting to c Rcpp if not sapply method

vector_var_matrix_calc_COR3 <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }

  p <- dim(MATR)[1]
  m <- p*(p-1)/2
  tocomp <- unlist(sapply(1:m, function(i) (i - 1)*m + i:m))

  real.cov2 <- function(q, MATR, p, m, cumsum) {
    t1 <- ceiling(q/m)
    t2 <- q %% m
    t2 <- m*(t2 == 0) + t2*(t2 != 0)
    i <- sum(cumsum < t1)
    j <- i + t1 - cumsum[i]
    k <- sum(cumsum < t2)
    l <- k + t2 - cumsum[k]

    MATRij <- MATR[i,j]
    MATRkl <- MATR[k,l]
    MATRik <- MATR[i,k]
    MATRil <- MATR[i,l]
    MATRjk <- MATR[j,k]
    MATRjl <- MATR[j,l]
    
    (MATRij*MATRkl/2) * (MATRik^2 + MATRil^2 + MATRjk^2 + MATRjl^2) -
      MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
      MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
      (MATRik*MATRjl + MATRil*MATRjk)
  }
  
  cumsum <- c(0, cumsum((p - 1):1))
  pelet <- mclapply(tocomp, real.cov2, MATR = MATR,
                   p = p, m = m, cumsum = cumsum,
                   mc.cores = ifelse(.Platform$OS.type == "windows", 1, ncores))
  pelet <- vector2triangle(unlist(pelet), diag = T)
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

profvis({
  tt1 <- Sys.time()
  pelet <- vector_var_matrix_calc_COR(MATR)
  tt1 <- Sys.time() - tt1
  tt3 <- Sys.time()
  pelet3 <- vector_var_matrix_calc_COR3(MATR)
  tt3 <- Sys.time() - tt3
  
})

identical(round(pelet, 2), round(pelet3 ,2))
tt1
tt3
