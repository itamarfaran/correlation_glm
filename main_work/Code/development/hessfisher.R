CovObj <- Pelet_Cov

hessian(x = CovObj$alpha,
        func = function(A) minusloglik(theta = CovObj$theta,
                                       alpha = A, linkFun = linkFun$FUN,
                                       sick.data = sickDat,
                                       effective.N = CovObj$Est_N))


tmp <- function(A, alpha, i){
  A <- replace(alpha, i, A)
  minusloglik(theta = CovObj$theta,
              alpha = A, linkFun = linkFun$FUN,
              sick.data = sickDat,
              effective.N = CovObj$Est_N)
}

i <- 1
tmp2 <- function(alpha, i) grad(tmp, alpha[i], alpha = alpha, i = i)
res <- grad(tmp2, CovObj$alpha, i = i)


fisherMatrHess[,1] / res
