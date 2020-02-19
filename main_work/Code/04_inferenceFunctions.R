compute_mu_alpha_jacobian <- function(theta, alpha, d = 1, linkFun){
  return(
    jacobian(
      func = function(A) triangle2vector(
        linkFun$FUN(
          t = theta,
          a = A,
          d = d
          )
        ),
      x = alpha
    )
  )
}

compute_gee_variance <- function(CovObj, sick.data, linkFun = linkFunctions$multiplicative_identity,
                                 dim_alpha = 1, reg_lambda = 0, reg_p = 2, est_mu = TRUE, ncores = 1){
  
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick.data))/2
  d <- length(CovObj$alpha)/p
  
  mu_alpha_jacobian <- compute_mu_alpha_jacobian(CovObj$theta, CovObj$alpha, d = d, linkFun) 

  g11 <- if(est_mu){
    triangle2vector(
      linkFun$FUN(
        t = CovObj$theta,
        a = CovObj$alpha,
        d = length(CovObj$alpha)/p
      )
    )
  } else {
    colMeans(sick.data)
  }
  
  residuals <- sick.data - rep(1, nrow(sick.data)) %o% g11
  cov_mat <- t(residuals) %*% residuals / (nrow(sick.data) - 1)
  
  Sigma <- vector_var_matrix_calc_COR_C(vector2triangle(colMeans(sick.data)))
  solve_Sigma <- solve(Sigma)
  
  I0 <- t(mu_alpha_jacobian) %*% solve_Sigma %*% mu_alpha_jacobian
  solve_I0 <- solve(I0)
  
  I1 <- t(mu_alpha_jacobian) %*% solve_Sigma %*% cov_mat %*% solve_Sigma %*% mu_alpha_jacobian
  
  res <- solve_I0 %*% I1 %*% solve_I0 / nrow(sick.data)
  
  return(res)
}

compute_fisher_by_grad <- function(CovObj, sick.data, linkFun = linkFunctions$multiplicative_identity,
                                   dim_alpha = 1, reg_lambda = 0, reg_p = 2, ncores = 1){
  
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick.data))/2
  d <- length(CovObj$alpha)/p
  
  mu_alpha_jacobian <- compute_mu_alpha_jacobian(CovObj$theta, CovObj$alpha, d = d, linkFun) 
  
  g11 <- triangle2vector(
    linkFun$FUN(
      t = CovObj$theta,
      a = CovObj$alpha,
      d = length(CovObj$alpha)/p
    )
  )
  
  residuals <- sick.data - rep(1, nrow(sick.data)) %o% g11
  
  Sigma <- vector_var_matrix_calc_COR_C(vector2triangle(colMeans(sick.data)))/CovObj$Est_N
  temp_equation <- t(mu_alpha_jacobian) %*% solve(Sigma) %*% t(residuals)
  res <- temp_equation %*% t(temp_equation)
  return(res)
}

compute_fisher_by_hess <- function(CovObj, sick.data, linkFun = linkFunctions$multiplicative_identity,
                                   dim_alpha = 1, reg_lambda = 0, reg_p = 2, ncores = 1){
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  Sigma <- vector_var_matrix_calc_COR_C(vector2triangle(colMeans(sick.data)))/CovObj$Est_N
  inv_sigma <- solve(Sigma)
  to_deriv <- function(alpha){
    sum_of_squares(
      alpha = alpha,
      theta = CovObj$theta,
      sick.data = sick.data,
      inv_sigma = inv_sigma,
      linkFun = linkFun)
  }
  fisher_mat <- hessian(to_deriv, CovObj$alpha)
  return(fisher_mat)
}

compute_sandwhich_fisher_variance <- function(CovObj, sick.data, linkFun = linkFunctions$multiplicative_identity,
                                              dim_alpha = 1, reg_lambda = 0, reg_p = 2, ncores = 1){
  grad_fisher <- compute_fisher_by_grad(
    CovObj = CovObj,
    sick.data = sick.data,
    linkFun = linkFun,
    dim_alpha = dim_alpha,
    reg_lambda = reg_lambda,
    reg_p = reg_p,
    ncores = ncores)
  hess_fisher <- compute_fisher_by_hess(
    CovObj = CovObj,
    sick.data = sick.data,
    linkFun = linkFun,
    dim_alpha = dim_alpha,
    reg_lambda = reg_lambda,
    reg_p = reg_p,
    ncores = ncores)
  hess_fisher_solve <- solve(hess_fisher)
  out <- hess_fisher_solve %*% grad_fisher %*% hess_fisher_solve
  return(out)
}

### todo: modify this according to hypothesis
### todo: I stopped reviewing here
build_hyp_test <- function(CovObj, VarMat, effectiveN, linkFun = linkFunctions$multiplicative_identity,
                           test = c("lower", "upper", "two-sided"),
                           sig.level = 0.05, p.adjust.method = p.adjust.methods, const = 1, Real){
  
  if(length(p.adjust.method) > 1) p.adjust.method <- p.adjust.method[1]
  if(length(test) > 1) test <- test[3]

  alpha_var_mat <- VarMat
  alpha_sd <- sqrt(diag(alpha_var_mat))
  
  dist_fun <- function(q) pnorm(q)
  critical_value <- qnorm(1 - sig.level/2)
  # critical_value <- qmvnorm(1 - sig.level/2, corr = cov2cor(alpha_var_mat))$quantile

  res <- data.frame(Est. = CovObj$alpha)
  res$Std. <- const*alpha_sd
  res$'Z-val' <- (res$Est. - linkFun$NULL_VAL)/res$Std.
  res$'Lower' <- res$Est. - critical_value*alpha_sd
  res$'Upper' <- res$Est. + critical_value*alpha_sd
  
  if(test == "lower"){ res$'P-val' <- round(dist_fun(res$'Z-val'), 5)
  }else if(test == "upper"){ res$'P-val' <- round(1 - dist_fun(res$'Z-val'), 5)
  }else{res$'P-val' <- round(2*(1-dist_fun(abs(res$'Z-val'))), 5)}
  
  res$'Adj.P-val' <- p.adjust(res$'P-val', method = p.adjust.method)
  res$Reject_H0 <- res$'Adj.P-val' < sig.level
  
  if(!missing(Real)) res$Real <- linkFun$INV(as.vector(Real))

  return(list(
    Results = res, Test = test, Significance = sig.level,
    p.adjust.method = p.adjust.method, CriticalVal = critical_value, Var_Mat = alpha_var_mat))
}

wilksTest <- function(covObj, healthy.dat, sick.dat, dim_alpha = 1, linkFun){
  if(class(healthy.dat) == "array") healthy.dat <- cor.matrix_to_norm.matrix(healthy.dat)
  if(class(sick.dat) == "array") sick.dat <- cor.matrix_to_norm.matrix(sick.dat)

  #Do a wilks test (chi-square)
  chisq <- -2*( minusloglik(theta = covObj$theta,
                            alpha = covObj$alpha,
                            linkFun = linkFun,
                            healthy.data = healthy.dat,
                            sick.data = sick.dat,
                            dim_alpha = dim_alpha) -
                  minusloglik(theta = colMeans(rbind(healthy.dat, sick.dat)),
                              alpha = rep(linkFun$NULL_VAL, length(covObj$alpha)),
                              linkFun = linkFun,
                              healthy.data = healthy.dat,
                              sick.data = sick.dat,
                              dim_alpha = dim_alpha) )
  
  c("Chisq_val" = chisq, "DF" = length(covObj$alpha),
    "Pval" = 1 - pchisq(chisq, length(covObj$alpha)))
}

multipleComparison <- function(healthy.data, sick.data,
                               p.adjust.method = p.adjust.methods, test = c("lower", "upper", "both")){
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  if(length(p.adjust.method > 1)) p.adjust.method <- p.adjust.method[1]
  if(length(test) > 1) test <- test[3]
  
  healthy_means <- apply(healthy.data, 2, mean)
  healthy_vars <- apply(healthy.data, 2, var)
  healthy_ns <- apply(healthy.data, 2, length)
  sick_means <- apply(sick.data, 2, mean)
  sick_vars <- apply(sick.data, 2, var)
  sick_ns <- apply(sick.data, 2, length)
  dfs <- (healthy_vars/healthy_ns + sick_vars/sick_ns)^2/
    ((healthy_vars/healthy_ns)^2/(healthy_ns - 1) + (sick_vars/sick_ns)^2/(sick_ns - 1))
  tvals <- (sick_means - healthy_means)/sqrt(healthy_vars/healthy_ns + sick_vars/sick_ns)
  
  if(test == "lower") Pvals <- pt(tvals, dfs)
  if(test == "upper") Pvals <- 1 - pt(tvals, dfs)
  if(test == "both") Pvals <- 2*pt(abs(tvals), dfs, lower.tail = F)
  p.adjust(Pvals, p.adjust.method)
}
