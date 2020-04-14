theta <- Pelet_Cov$theta
alpha <- Pelet_Cov$alpha
eff_N <- Pelet_Cov$Est_N
# linkFun <- linkFunctions$Benjamini$FUN
linkFun <- linkFunctions$Benjamini

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
}

loglikgrad_uni <- function(obs, CovObj, linkFun){
  x <- grad(function(x) loglik_uni(obs = obs, theta = CovObj$theta, alpha = x, Eff.N = CovObj$Est_N, linkFun = linkFun),
            CovObj$alpha)
  return(x %*% t(x))
}

computeBmatr <- function(CovObj, sickDat, silent = FALSE, linkFun, ncores = 1){
  Bmatr <- mclapply(1:nrow(sickDat), function(j) loglikgrad_uni(sickDat[j,], CovObj, linkFun), mc.cores = ncores)
  Bmatr <- Bmatr %>% (function(list){
    L <- length(list)
    pelet <- matrix(0, nrow = nrow(list[[1]]), ncol = ncol(list[[1]]))
    for(i in 1:L) pelet <- pelet + list[[i]]
    pelet
  })
  return(Bmatr)
}

computeFisherByGrad <- function(CovObj, sickDat, linkFun, U1 = TRUE, ncores = 1){
  gamma <- function(theta, alpha, eff_N, linkFun){
    forGrad <- function(A, i){
      A <- replace(alpha, i, A)
      g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(A))))
      g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11))
      e21 <- eigen(g21/eff_N, symmetric = TRUE)
      U1 <- e21$vectors
      D1 <- e21$values
      
      return(sum(log(D1)) + t(g11) %*% U1 %*% diag(1/D1) %*% t(U1) %*% g11)
    }
    mclapply(1:length(alpha), function(i) grad(func = forGrad, x = alpha[i], i = i), mc.cores = ncores) %>% unlist()
  }
  kappa <- function(x, theta, alpha, eff_N, linkFun, U1 = TRUE){
    if(U1){
      U1 <- NULL
    } else {
      e21 <- eigen(g21/n.effective_D, symmetric = TRUE)
      U1 <- e21$vectors
    }
    forGrad <- function(A){
      g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(A))))
      g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11))
      if(is.null(U1)){
        SigSolve <- solve(g21)
      } else {
        e21 <- eigen(g21/eff_N, symmetric = TRUE, only.values = T)
        D1 <- e21$values
        
        SigSolve <- U1 %*% diag(1/D1) %*% t(U1)
      }
      return(SigSolve %*% (x - 2*g11) )
    }
    t(x) %*% jacobian(forGrad, alpha)
  }
  
  gammares <- gamma(theta = CovObj$theta, alpha = CovObj$alpha, eff_N = CovObj$Est_N, linkFun = linkFun)
  kappares <- mclapply(1:nrow(sickDat),
                       function(i) kappa(x = sickDat[i,], theta = CovObj$theta, alpha = CovObj$alpha, eff_N = CovObj$Est_N,
                                         linkFun = linkFun, U1 = U1),
                       mc.cores = ncores)
  
  kappares <- do.call(rbind, kappares)
  kapgam <- gammares %o% colSums(kappares)
  kapkap <- mclapply(1:nrow(sickDat), function(i) kappares[i,] %o% kappares[i,], mc.cores = ncores) %>%
    simplify2array() %>% calculate_mean_matrix(do.mean = FALSE)
  
  return(0.25*(nrow(sickDat) * gammares %o% gammares + kapgam + t(kapgam) + kapkap))
}

# profvis({
tt <- Sys.time()
resnew <- computeFisherByGrad(Pelet_Cov, cor.matrix_to_norm.matrix(sampleData$sick), linkFun = linkFun$FUN, U1 = T, ncores = ncores)
tt <- Sys.time() - tt

tt2 <- Sys.time()
resold <- computeBmatr(Pelet_Cov, cor.matrix_to_norm.matrix(sampleData$sick), linkFun = linkFun$FUN, ncores = ncores)
tt2 <- Sys.time() - tt2

# })
