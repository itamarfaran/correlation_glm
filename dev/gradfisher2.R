library(microbenchmark)

linkFun <- linkFunctions$Benjamini

computeFisherByGrad <- function(CovObj, sickDat, linkFun, U1 = TRUE, ncores = 1){
  U1 = TRUE
  # U1 == TRUE    =>   Compute U1 every differntiation
  # U1 == FALSE   =>   Compute U1 once, before differntiation
  # U1 == matrix  =>   use given U1
  
  if(is.logical(U1)){
    if(U1) {
      U1 = NULL
    } else {
      g11 <- as.matrix(CovObj$theta*triangle2vector(create_alpha_mat(linkFun(CovObj$alpha))))
      g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11))/CovObj$Est_N
      e21 <- eigen(g21, symmetric = TRUE)
      U1 <- e21$vectors
    }
  }

  gamma <- function(theta, alpha, eff_N, linkFun, U1){
    forGrad <- function(A, i){
      A <- replace(alpha, i, A)
      g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(A))))
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
    mclapply(seq_along(alpha), function(i) grad(func = forGrad, x = alpha[i], i = i), mc.cores = ncores) %>% unlist()
  }
  kappa <- function(x, theta, alpha, eff_N, linkFun, U1){
    forGrad <- function(A){
      g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(A))))
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
  
  gammares <- gamma(theta = CovObj$theta, alpha = CovObj$alpha, eff_N = CovObj$Est_N, linkFun = linkFun, U1 = U1)
  kappares <- do.call(rbind, mclapply(
    seq_len(nrow(sickDat)),
    function(i) kappa(x = sickDat[i,], theta = CovObj$theta, alpha = CovObj$alpha, eff_N = CovObj$Est_N,
                      linkFun = linkFun, U1 = U1),
    mc.cores = ncores))
  
  kapgam <- gammares %o% colSums(kappares)
  kapkap <- mclapply(seq_len(nrow(sickDat)), function(i) kappares[i,] %o% kappares[i,], mc.cores = ncores) %>%
    simplify2array() %>% calculate_mean_matrix(do.mean = FALSE)
  
  return(0.25*(nrow(sickDat) * gammares %o% gammares + kapgam + t(kapgam) + kapkap))
}

loglik_uni <- function(obs, theta, alpha, Eff.N, linkFun){
  
  if(missing(alpha)){
    meanMat <- vector2triangle(theta)
  } else {
    meanMat <- vector2triangle(theta)*create_alpha_mat(linkFun(alpha))
  }
  varMat <- vector_var_matrix_calc_COR_C(meanMat)/Eff.N
  meanVect <- triangle2vector(meanMat)
  
  eigenDec <- eigen(varMat, symmetric = TRUE)
  U <- eigenDec$vectors
  D <- eigenDec$values
  dists <- t(U) %*% (obs - meanVect)
  -0.5*(sum(log(D)) + sum(dists^2/D))
  # new_obs <- t(U) %*% obs
  # new_meanVect <- t(U) %*% meanVect
  # new_invD <- diag(1/D)
  # -0.5*(t(new_obs) %*% new_invD %*% new_obs - 2*t(new_meanVect) %*% new_invD %*% new_obs)
}

loglikgrad_uni <- function(obs, CovObj, linkFun){
  x <- grad(function(x) loglik_uni(obs = obs, theta = CovObj$theta, alpha = x, Eff.N = CovObj$Est_N, linkFun = linkFun),
            CovObj$alpha)
  return(x %*% t(x))
}

computeBmatr <- function(CovObj, sickDat, silent = FALSE, linkFun, ncores = 1){
  Bmatr <- mclapply(seq_len(nrow(sickDat)), function(j) loglikgrad_uni(sickDat[j,], CovObj, linkFun), mc.cores = ncores)
  Bmatr <- Bmatr %>% (function(list){
    L <- length(list)
    pelet <- matrix(0, nrow = nrow(list[[1]]), ncol = ncol(list[[1]]))
    for(i in 1:L) pelet <- pelet + list[[i]]
    pelet
  })
  return(Bmatr)
}


# res_mbm =
#   microbenchmark(
#     computeBmatr(Pelet_Cov, cor.matrix_to_norm.matrix(sampleData$sick), linkFun = linkFun$FUN),
#     computeFisherByGrad(Pelet_Cov, cor.matrix_to_norm.matrix(sampleData$sick), linkFun = linkFun$FUN),
#     times = 10
#   )
# 
# res_mbm
# res_mbm =
#   microbenchmark(
#     computeFisherByGrad(Pelet_Cov, cor.matrix_to_norm.matrix(sampleData$sick), linkFun = linkFun$FUN, U1 = TRUE),
#     computeFisherByGrad(Pelet_Cov, cor.matrix_to_norm.matrix(sampleData$sick), linkFun = linkFun$FUN, U1 = FALSE),
#     times = 10
#   )
# 
# res_mbm

# CovObj = Pelet_Cov
sickDat = cor.matrix_to_norm.matrix(sampleData$sick)
# linkFun = linkFun$FUN

computeBmatr(Pelet_Cov, sickDat, linkFun = linkFun$FUN)
computeFisherByGrad(Pelet_Cov, sickDat, linkFun = linkFun$FUN)

