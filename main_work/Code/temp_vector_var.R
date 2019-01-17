source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions2.R")

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

p <- 10
MATR <- build_parameters(p, 0.5, c(0,1))$Corr.mat
# MATR <- matrix(1:9, ncol = 3)
# MATR <- MATR + t(MATR) + diag(3)*9

vector_var_matrix_calc_COR <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i)))
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p))
  
  pelet <- matrix(0, nrow = m, ncol = m)
  for(i1 in 1:m){
    for(j1 in i1:m){
      i <- order_vecti[i1]
      j <- order_vectj[i1]
      k <- order_vecti[j1]
      l <- order_vectj[j1]
      
      MATRij <- MATR[i,j]
      MATRkl <- MATR[k,l]
      MATRik <- MATR[i,k]
      MATRil <- MATR[i,l]
      MATRjk <- MATR[j,k]
      MATRjl <- MATR[j,l]
      
      pelet[i1,j1] <-
        (MATRij*MATRkl/2) * (MATRik^2 + MATRil^2 + MATRjk^2 + MATRjl^2) -
        MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
        MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
        (MATRik*MATRjl + MATRil*MATRjk)
    }
  }
  
  pelet <- pelet + t(pelet) - diag(diag(pelet))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

cppFunction(
'NumericMatrix corcalc_c(NumericMatrix MATR,
  int p, int m, NumericVector order_vecti, NumericVector order_vectj) {
  
  NumericMatrix pelet(m, m); 

  for (int i1 = 0; i1 < m; i1++) {
    for (int j1 = 0; j1 < m; j1++) {
      int i = order_vecti[i1];
      int j = order_vectj[i1];
      int k = order_vecti[j1];
      int l = order_vectj[j1];

      int MATRij = MATR(i,j);
      int MATRkl = MATR(k,l);
      int MATRik = MATR(i,k);
      int MATRil = MATR(i,l);
      int MATRjk = MATR(j,k);
      int MATRjl = MATR(j,l);

      pelet(i1,j1) =
        (MATRij*MATRkl/2) * (pow(MATRik, 2) + pow(MATRil, 2) + pow(MATRjk, 2) + pow(MATRjl, 2)) -
        MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
        MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
        (MATRik*MATRjl + MATRil*MATRjk);
    }
  }
  return(pelet);
}')

corcalc_R <- function(MATR, p, m, order_vecti, order_vectj){
  pelet <- matrix(0, nrow = m, ncol = m)
  for(i1 in 1:m){
    for(j1 in i1:m){
      i <- order_vecti[i1]
      j <- order_vectj[i1]
      k <- order_vecti[j1]
      l <- order_vectj[j1]
      
      MATRij <- MATR[i,j]
      MATRkl <- MATR[k,l]
      MATRik <- MATR[i,k]
      MATRil <- MATR[i,l]
      MATRjk <- MATR[j,k]
      MATRjl <- MATR[j,l]
      
      pelet[i1,j1] <-
        (MATRij*MATRkl/2) * (MATRik^2 + MATRil^2 + MATRjk^2 + MATRjl^2) -
        MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
        MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
        (MATRik*MATRjl + MATRil*MATRjk)
    }
  }
  return(pelet)
}

vector_var_matrix_calc_COR_CR <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                          reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i)))
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p))
  
  pelet <- corcalc_R(MATR, p, m, order_vecti, order_vectj)
  pelet <- pelet + t(pelet) - diag(diag(pelet))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
} # todo : try exporting to c Rcpp if not sapply method

vector_var_matrix_calc_COR_C <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                          reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i))) - 1
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p)) - 1
  
  pelet <- corcalc_c(MATR, p, m, order_vecti, order_vectj)
  pelet <- pelet + t(pelet) - diag(diag(pelet))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
} # todo : try exporting to c Rcpp if not sapply method

vector_var_matrix_calc_COR_par <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
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
  if(reg_par != 0) pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

# profvis({
  tt1 <- Sys.time()
  pelet1 <- vector_var_matrix_calc_COR(MATR)
  tt1 <- Sys.time() - tt1
  tt2 <- Sys.time()
  pelet2 <- vector_var_matrix_calc_COR_C(MATR)
  tt2 <- Sys.time() - tt2
  tt3 <- Sys.time()
  #pelet3 <- vector_var_matrix_calc_COR_par(MATR)
  tt3 <- Sys.time() - tt3
# })

identical(round(pelet1, 2), round(pelet2 ,2))
#identical(round(pelet2, 2), round(pelet3 ,2))
tt1
tt2
tt3
#tt4

# rm(pelet1, pelet2, pelet3, pelet4)
gc()
