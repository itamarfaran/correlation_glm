CovObj <- Pelet_Cov
sickDat <- cor.matrix_to_norm.matrix(sampleData$sick)

loglik_uni <- function(obs, theta, alpha, Eff.N){
  
  if(missing(alpha)){
    meanMat <- vector2triangle(theta)
  } else {
    meanMat <- vector2triangle(theta)*create_alpha_mat(alpha)
  }
  varMat <- vector_var_matrix_calc_COR_C(meanMat)/Eff.N
  meanVect <- triangle2vector(meanMat)
  
  eigenDec <- eigen(varMat, symmetric = TRUE)
  U <- eigenDec$vectors
  D <- eigenDec$values
  
  dists <- t(U) %*% (obs - meanVect)
  
  -0.5*(sum(log(D)) + sum(dists^2/D))
}

loglikgrad_uni <- function(obs, CovObj){
  grad(function(x) loglik_uni(obs = obs, theta = CovObj$theta, alpha = x, Eff.N = CovObj$Est_N), CovObj$alpha) %>%
    (function(x) x %*% t(x))
}

computeBmatr <- function(CovObj, sickDat, silent = FALSE, ncores = .GlobalEnv$ncores){
  if(ncores == 1) {
    Bmatr <- lapply(1:nrow(sickDat), function(j) loglikgrad_uni(sickDat[j,], CovObj))
  } else {
    rawFun <- function(j) loglikgrad_uni(sickDat[j,], CovObj)
    cl <<- makeCluster(ncores)
    if(!silent) message("In 'computeBmatr': Cluster 'cl' opened, saved to global environment.")
    clusterEvalQ(cl, library(mvtnorm))
    clusterEvalQ(cl, library(dplyr))
    clusterEvalQ(cl, library(numDeriv))
    clusterEvalQ(cl, library(matrixcalc))
    clusterExport(cl, c("sickDat", "CovObj"), envir = environment())
    clusterExport(cl, c("loglik_uni", "loglikgrad_uni", "vector2triangle","create_alpha_mat",
                        "vector_var_matrix_calc_COR_C", "triangle2vector", "corcalc_c"), envir = .GlobalEnv)
    Bmatr <- parLapply(cl = cl, 1:nrow(sickDat), rawFun)
    terminateCL(silent)
  }
  Bmatr <- Bmatr %>% (function(list){
    L <- length(list)
    pelet <- matrix(0, nrow = nrow(list[[1]]), ncol = ncol(list[[1]]))
    for(i in 1:L) pelet <- pelet + list[[i]]
    pelet
  })
  
  return(Bmatr)
}

N <- nrow(sickDat)
p <- length(CovObj$alpha)
m <- 0.5*p*(p-1)

logdet <- function(theta, alpha, Eff.N){
  grad(function(x)
    sum(log(eigen(
        vector_var_matrix_calc_COR_C(
          vector2triangle(theta)*create_alpha_mat(x)/Eff.N),
        symmetric = T, only.values = T
      )$values)), alpha)
}
psi <- logdet(CovObj$theta, CovObj$alpha, CovObj$Est_N)

xi_fun <- function(theta, alpha, Eff.N, sickDati){
  dist <- sickDati - triangle2vector(vector2triangle(theta)*create_alpha_mat(alpha))
  g2s <- solve(vector_var_matrix_calc_COR_C(vector2triangle(theta)*create_alpha_mat(alpha)/Eff.N))
  t(dist) %*% g2s %*% dist
}
xi <- sapply(1:nrow(sickDat),
             function(i)
               grad(function(x) xi_fun(CovObj$theta, x, CovObj$Est_N, sickDat[i,]),
                    CovObj$alpha))

sumxi <- rowSums(xi)
psipsi <- N * (psi %o% psi)
xipsi <- sumxi %o% psi
tmp <- sapply(1:ncol(xi), function(i) xi[,i] %o% xi[,i], simplify = "array")
xixi <- matrix(0, p, p)
for(i in 1:N) xixi <- xixi + tmp[,,i]

(psipsi + xipsi + t(xipsi) + xixi)*0.25/res

res <- computeBmatr(CovObj, sickDat)

dim(xi)
g2 <- vector_var_matrix_calc_COR_C(vector2triangle(CovObj$theta)*create_alpha_mat(CovObj$alpha)/CovObj$Est_N)
g2slv <- solve(g2)

tmp <- lapply(1:(m^2), function(i)
  grad(function(x)
    vector_var_matrix_calc_COR_C(
      vector2triangle(CovObj$theta)*create_alpha_mat(x)/CovObj$Est_N)[i], CovObj$alpha))
dg2 <- array(dim = c(m, m, p))
for(i in 1:p) dg2[,,i] <- tmp[[i]]

lapply(1:p, function(i) -g2slv*dg2[,,i]*g2slv)


length(vector2triangle(CovObj$theta)*create_alpha_mat(CovObj$alpha)/CovObj$Est_N)

length(sol)


dim(ans)
ans[,,1]
