loglik_uni <- function(obs, theta, alpha = NULL, Eff.N){
  
  if(length(alpha) == 0){
    meanMat <- vector_to_triangle(theta)
  } else {
    meanMat <- vector_to_triangle(theta)*create_alpha_mat(alpha)
  }
  varMat <- vector_var_matrix_calc_COR(meanMat)/Eff.N
  meanVect <- triangle_to_vector(meanMat)
  
  eigenDec <- eigen(varMat)
  U <- eigenDec$vectors
  D <- eigenDec$values
  
  dists <- t(U) %*% (obs - meanVect)
  
  -0.5*(sum(log(D)) + sum(dists^2/D))
}

loglikgrad_uni <- function(obs, CovObj){
  grad(function(x) loglik_uni(obs = obs, theta = CovObj$theta, alpha = x, Eff.N = CovObj$Est_N), CovObj$alpha) %>%
    (function(x) x %*% t(x))
}

computeBmatr <- function(sickDat, CovObj){
  Bmatr <- lapply(1:dim(sickDat)[3], function(j) loglikgrad_uni(cor.matrix_to_norm.matrix(sickDat)[j,], CovObj)) %>%
    (function(list){
      L <- length(list)
      pelet <- matrix(0, nrow = nrow(list[[1]]), ncol = ncol(list[[1]]))
      for(i in 1:L) pelet <- pelet + list[[i]]
      pelet
    }) #/dim(sickDat)[3]
  
  return(Bmatr)
}

# Amatr <- b1Hess/dim(sick_tmp)[3]
# 
# Amatr
# Bmatrix
# 
# plot(Amatr, Bmatrix)
# 
# Cmat <- solve(Amatr) %>% (function(X) X %*% Bmatr %*% X) 
#  
# sqrt(diag(Cmat))/sqrt(diag(solve(Amatr)))
# 
# is.positive.definite(force_symmetry(solve(Amatr) - Cmat))