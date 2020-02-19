computeFisherByGrad_old <- function(CovObj, sickDat, linkFun = linkFunctions$multiplicative_identity,
                                    dim_alpha = 1, nonpositive = "Stop", reg_lambda = 0, U1 = TRUE, ncores = 1){
  U1 = TRUE
  # U1 == TRUE    =>   Compute U1 every diffrentiation
  # U1 == FALSE   =>   Compute U1 once, before diffrentiation
  # U1 == matrix  =>   use given U1
  
  if(is.logical(U1)){
    if(U1) {
      U1 = NULL
    } else {
      g11 <- as.matrix(triangle2vector(linkFun$FUN(t = CovObj$theta, a = CovObj$alpha, d = dim_alpha))) 
      g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11, nonpositive = nonpositive))/CovObj$Est_N
      e21 <- eigen(g21, symmetric = TRUE)
      U1 <- e21$vectors
    }
  }
  
  gamma <- function(theta, alpha, eff_N, linkFun, dim_alpha, U1){
    forGrad <- function(A, i){
      A <- replace(alpha, i, A)
      g11 <- as.matrix(triangle2vector(linkFun$FUN(t = theta, a = A, d = dim_alpha))) 
      g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11))/eff_N
      if(is.null(U1)){
        e21 <- eigen(g21, symmetric = TRUE)
        U1 <- e21$vectors
      } else {
        e21 <- eigen(g21, symmetric = TRUE, only.values = TRUE)
      }
      D1 <- e21$values
      
      return(sum(log(D1)) + t(g11) %*% U1 %*% diag(1/D1) %*% t(U1) %*% g11)
    }
    unlist(mclapply(
      1:length(alpha), function(i) grad(func = forGrad, x = alpha[i], i = i), mc.cores = ncores
    ))
  }
  kappa <- function(x, theta, alpha, eff_N, linkFun, dim_alpha, U1){
    forGrad <- function(A){
      g11 <- as.matrix(triangle2vector(linkFun$FUN(t = theta, a = A, d = dim_alpha))) 
      g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11))/eff_N
      if(is.null(U1)){
        sig_solve <- solve(g21)
      } else {
        D1 <- eigen(g21, symmetric = TRUE, only.values = TRUE)$values
        sig_solve <- U1 %*% diag(1/D1) %*% t(U1)
      }
      return(t(x) %*% sig_solve %*% (x - 2*g11) )
    }
    grad(forGrad, alpha)
  }
  
  gammares <- gamma(theta = CovObj$theta, alpha = CovObj$alpha, eff_N = CovObj$Est_N, linkFun = linkFun, dim_alpha = dim_alpha, U1 = U1)
  kappares <- do.call(rbind, mclapply(
    1:nrow(sickDat),
    function(i) kappa(x = sickDat[i,], theta = CovObj$theta, alpha = CovObj$alpha, eff_N = CovObj$Est_N,
                      linkFun = linkFun, dim_alpha = dim_alpha, U1 = U1),
    mc.cores = ncores))
  
  kapgam <- gammares %o% colSums(kappares)
  kapkap <- mclapply(1:nrow(sickDat), function(i) kappares[i,] %o% kappares[i,], mc.cores = ncores) %>%
    simplify2array() %>% calculate_mean_matrix(do.mean = FALSE)
  
  return(0.25*(nrow(sickDat) * gammares %o% gammares + kapgam + t(kapgam) + kapkap) + 2*reg_lambda*CovObj$alpha)
}

computeFisherByHess_old <- function(CovObj, sickDat, method = c("Hess", "Grad"), linkFun, dim_alpha = 1,
                                    nonpositive = "Stop", reg_lambda = 0, ncores = 1, silent = FALSE){
  if(class(sickDat) == "array") sickDat <- cor.matrix_to_norm.matrix(sickDat) 
  
  method <- method[1]
  output <- # diag(rep(2*reg_lambda, length(CovObj$alpha))) +
    hessian(
      x = CovObj$alpha,
      func = function(A) minusloglik(theta = CovObj$theta,
                                     alpha = A, linkFun = linkFun,
                                     sick.data = sickDat,
                                     effective.N = CovObj$Est_N,
                                     dim_alpha = dim_alpha,
                                     nonpositive = nonpositive))
  
  return(output)
}

