computeFisherByGrad <- function(CovObj, sickDat, linkFun, dim_alpha = 1, U1 = TRUE, ncores = 1){
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
  
  gamma <- function(theta, alpha, eff_N, linkFun, dim_alpha, U1){
    forGrad <- function(A, i){
      A <- replace(alpha, i, A)
      g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(A), dim_alpha = dim_alpha)))
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
      g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(A), dim_alpha = dim_alpha)))
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
  
  return(0.25*(nrow(sickDat) * gammares %o% gammares + kapgam + t(kapgam) + kapkap))
}

ComputeFisher <- function(CovObj, sickDat, method = c("Hess", "Grad"), linkFun, dim_alpha = 1, ncores = 1, silent = FALSE){
  if(class(sickDat) == "array") sickDat <- cor.matrix_to_norm.matrix(sickDat) 
  if(missing(linkFun)) linkFun <- list(FUN = function(x) x, INV = function(x) x)
  
  method <- method[1]
  if(method == "Hess") output <- hessian(x = CovObj$alpha,
                                        func = function(A) minusloglik(theta = CovObj$theta,
                                                                       alpha = A, linkFun = linkFun$FUN,
                                                                       sick.data = sickDat,
                                                                       effective.N = CovObj$Est_N, 
                                                                       dim_alpha = dim_alpha))
  if(method == "Grad") output <- computeFisherByGrad(CovObj, sickDat, linkFun = linkFun$FUN, dim_alpha = dim_alpha, ncores = ncores)
  
  return(output)
}

### todo: modify this according to hypothesis
build_hyp.test <- function(CovObj, FisherMatr, effectiveN, linkFun, test = c("lower", "upper", "two-sided"),
                             sig.level = 0.05, p.adjust.method = p.adjust.methods, const = 1, Real){
  
  if(length(p.adjust.method) > 1) p.adjust.method <- p.adjust.method[1]
  if(length(test) > 1) test <- test[3]
  if(missing(linkFun)) linkFun <- list(FUN = function(x) x, INV = function(x) x)

  alpha_var_mat <- solve(FisherMatr)
  alpha_sd <- sqrt(diag(alpha_var_mat))
  
  # if(!missing(effectiveN)){
  #   dist_fun <- function(q) pt(q, effectiveN)
  #   critical_value <- qmvt(1 - sig.level/2, corr = cov2cor(alpha_var_mat), df = ceiling(effectiveN))$quantile
  # } else {
    dist_fun <- function(q) pnorm(q)
    critical_value <- qnorm(1 - sig.level/2)#qmvnorm(1 - sig.level/2, corr = cov2cor(alpha_var_mat))$quantile
  # }
  
  res <- data.frame(Est. = CovObj$alpha)
  res$Std. <- const*alpha_sd
  res$'Z-val' <- (res$Est. - linkFun$INV(1))/res$Std.

  res$'Lower' <- res$Est. - critical_value*alpha_sd
  res$'Upper' <- res$Est. + critical_value*alpha_sd
  
  if(test == "lower"){ res$'P-val' <- round(dist_fun(res$'Z-val'),5)
  }else if(test == "upper"){ res$'P-val' <- round(1 - dist_fun(res$'Z-val'),5)
  }else{res$'P-val' <- round(2*(1-dist_fun(abs(res$'Z-val'))),5)}
  
  res$'Adj.P-val' <- p.adjust(res$'P-val', method = p.adjust.method)
  res$Reject_H0 <- res$'Adj.P-val' < sig.level
  
  if(!missing(Real)) res$Real <- linkFun$INV(as.vector(Real))
  # if(!missing(effectiveN)){
  #   colnames(res)[colnames(res) == "Z-val"] <- "T-val"
  # }
  
  return(list(Results = res, Test = test, Significance = sig.level,
              p.adjust.method = p.adjust.method, CriticalVal = critical_value, Var_Mat = alpha_var_mat))
}

wilksTest <- function(covObj, healthy.dat, sick.dat, dim_alpha = 1, linkFun){
  if(class(healthy.dat) == "array") healthy.dat <- cor.matrix_to_norm.matrix(healthy.dat)
  if(class(sick.dat) == "array") sick.dat <- cor.matrix_to_norm.matrix(sick.dat)
  if(missing(linkFun)) linkFun <- list(FUN = function(x) x, INV = function(x) x)
  
  #Do a wilks test (chi-square)
  chisq <- -2*( minusloglik(theta = covObj$theta,
                            alpha = linkFun$FUN(covObj$alpha),
                            healthy.data = healthy.dat,
                            sick.data = sick.dat,
                            dim_alpha = dim_alpha) -
                  minusloglik(theta = rbind(healthy.dat, sick.dat) %>% colMeans(),
                              alpha = rep(1, length(covObj$alpha)),
                              healthy.data = healthy.dat,
                              sick.data = sick.dat,
                              dim_alpha = dim_alpha) )
  
  c("Chisq_val" = chisq, "DF" = length(covObj$alpha),
    "Pval" = 1 - pchisq(chisq, length(covObj$alpha)))
}
