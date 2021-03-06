loglik_uni <- function(obs, theta, alpha, Eff.N, linkFun){
  
  if(missing(alpha)){
    meanMat <- vector2triangle(theta)
  } else {
    meanMat <- vector2triangle(theta)*create_alpha_mat(linkFun(alpha))
    varMat <- vector_var_matrix_calc_COR_C(meanMat)/Eff.N
    meanVect <- triangle2vector(meanMat)
    
    eigenDec <- eigen(varMat, symmetric = TRUE)
    U <- eigenDec$vectors
    D <- eigenDec$values
    dists <- t(U) %*% (obs - meanVect)
    -0.5*(sum(log(D)) + sum(dists^2/D))
  }
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

ComputeFisher <- function(CovObj, sickDat, method = c("Hess", "Grad"), linkFun, ncores = 1, silent = FALSE){
  if(class(sickDat) == "array") sickDat <- cor.matrix_to_norm.matrix(sickDat) 
  if(missing(linkFun)) linkFun <- list(FUN = function(x) x, INV = function(x) x)
  
  method <- method[1]
  if(method == "Hess") pelet <- hessian(x = CovObj$alpha,
                                        func = function(A) minusloglik(theta = CovObj$theta,
                                                                       alpha = A, linkFun = linkFun$FUN,
                                                                       sick.data = sickDat,
                                                                       effective.N = CovObj$Est_N))
  if(method == "Grad") pelet <- computeBmatr(CovObj, sickDat, silent = silent, linkFun = linkFun$FUN, ncores = ncores)
  
  return(pelet)
}

build_hyp.test <- function(Estimate.Loop2_object, FisherMatr, effectiveN, linkFun, test = c("lower", "upper", "two-sided"),
                             sig.level = 0.05, p.adjust.method = p.adjust.methods, const = 1, Real){
  
  if(length(p.adjust.method) > 1) p.adjust.method <- p.adjust.method[1]
  if(length(test) > 1) test <- test[3]
  if(missing(linkFun)) linkFun <- list(FUN = function(x) x, INV = function(x) x)
  obj <- Estimate.Loop2_object
  
  alpha_var_mat <- solve(FisherMatr)
  alpha_sd <- sqrt(diag(alpha_var_mat))
  
  # if(!missing(effectiveN)){
  #   dist_fun <- function(q) pt(q, effectiveN)
  #   critical_value <- qmvt(1 - sig.level/2, corr = cov2cor(alpha_var_mat), df = ceiling(effectiveN))$quantile
  # } else {
    dist_fun <- function(q) pnorm(q)
    critical_value <- qmvnorm(1 - sig.level/2, corr = cov2cor(alpha_var_mat))$quantile
  # }
  
  res <- data.frame(Est. = obj$alpha)
  res$Std. <- const*alpha_sd
  res$'Z-val' <- (res$Est. - linkFun$INV(1))/res$Std.

  res$'Lower' <- res$Est. - critical_value*alpha_sd
  res$'Upper' <- res$Est. + critical_value*alpha_sd
  
  if(test == "lower"){ res$'P-val' <- round(dist_fun(res$'Z-val'),5)
  }else if(test == "upper"){ res$'P-val' <- round(1 - dist_fun(res$'Z-val'),5)
  }else{res$'P-val' <- round(2*(1-dist_fun(abs(res$'Z-val'))),5)}
  
  res$'Adj.P-val' <- p.adjust(res$'P-val', method = p.adjust.method)
  res$Reject_H0 <- res$'Adj.P-val' < sig.level
  
  if(!missing(Real)) res$Real <- linkFun$INV(Real)
  # if(!missing(effectiveN)){
  #   colnames(res)[colnames(res) == "Z-val"] <- "T-val"
  # }
  
  return(list(Results = res, Test = test, Significance = sig.level,
              p.adjust.method = p.adjust.method, CriticalVal = critical_value, Var_Mat = alpha_var_mat))
}

wilksTest <- function(covObj, healthy.dat, sick.dat, linkFun){
  if(class(healthy.dat) == "array") healthy.dat <- cor.matrix_to_norm.matrix(healthy.dat)
  if(class(sick.dat) == "array") sick.dat <- cor.matrix_to_norm.matrix(sick.dat)
  if(missing(linkFun)) linkFun <- list(FUN = function(x) x, INV = function(x) x)
  
  #Do a wilks test (chi-square)
  chisq <- -2*( minusloglik(theta = covObj$theta,
                            alpha = linkFun$FUN(covObj$alpha),
                            healthy.data = healthy.dat,
                            sick.data = sick.dat) -
                  minusloglik(theta = rbind(healthy.dat, sick.dat) %>% colMeans(),
                              alpha = rep(1, length(covObj$alpha)),
                              healthy.data = healthy.dat,
                              sick.data = sick.dat) )
  
  c("Chisq_val" = chisq, "DF" = length(covObj$alpha),
    "Pval" = 1 - pchisq(chisq, length(covObj$alpha)))
  
}