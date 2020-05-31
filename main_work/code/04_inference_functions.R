compute_mu_alpha_jacobian <- function(type, alpha, healthy_dt, sick_dt, d = 1, linkFun){
  func <- if(type == 'sick'){
    function(A) triangle2vector(
      linkFun$FUN(
        t = theta_of_alpha(A, healthy_dt, sick_dt, linkFun = linkFun, d = d),
        a = A,
        d = d
      )
    )
  } else if(type == 'healthy') {
    function(A) theta_of_alpha(A, healthy_dt, sick_dt, linkFun = linkFun, d = d)
  }
  return(
    jacobian(func = func, x = alpha)
  )
}

efrons_rms <- function(m, p = NULL){
  if(class(m) != 'matrix')
    m <- vector2triangle(m, diag_value = 1)
  if(!is.square.matrix(m))
    stop('m is not square')
  if(!is.positive.semi.definite(m))
    stop('m is not semi positive definite')
  if(!all(diag(m) == 1))
    stop('diag of m is not 1')
  
  m_vect <- triangle2vector(m, diag = FALSE)
  sum_sqrd_corrs <- sum(m_vect^2)
  rms <- sqrt(sum_sqrd_corrs/choose(ncol(m), 2))
  
  if(!is.null(p))
    rms <- sqrt(p/(p - 1) * (rms^2 - 1/(p - 1)))
  
  return(rms)
}

efrons_effective_sample_size <- function(n, rms){
  out <- n/(1 + (n - 1) * rms^2)
  return(out)
}

compute_gee_variance <- function(
  cov_obj, healthy_dt, sick_dt, est_mu = TRUE,
  reg_lambda = 0, reg_p = 2){
  
  # compute n/d factor <- function(){}
  
  compute_gee_raw <- function(type, list_){
    if(type == 'I0'){
      out <- t(list_$jacobian) %*% list_$solve_Sigma %*% list_$jacobian
    } else if (type == 'I1'){
      residuals <- list_$data - rep(1, nrow(list_$data)) %o% list_$expected_value
      cov_mat <- t(residuals) %*% residuals / list_$df # todo: add compute_estimated_n here
      out <- t(list_$jacobian) %*% list_$solve_Sigma %*% cov_mat %*% list_$solve_Sigma %*% list_$jacobian
    }
    out <- out*nrow(list_$data)
    return(out)
  }
  
  linkFun <- cov_obj$linkFun
  
  healthy_data <- convert_corr_array_to_data_matrix_test(healthy_dt)
  sick_data <- convert_corr_array_to_data_matrix_test(sick_dt)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick_data))/2
  d <- length(cov_obj$alpha)/p
  
  healthy_list <- list(
    data = healthy_data,
    jacobian = compute_mu_alpha_jacobian(
      type = 'healthy',
      alpha = cov_obj$alpha,
      healthy_dt = healthy_data,
      sick_dt = sick_data,
      d = d,
      linkFun = linkFun),
    expected_value = if(est_mu) cov_obj$theta else colMeans(healthy_data),
    solve_Sigma = solve(corrmat_covariance_from_dt(healthy_data, est_n = T)),
    df = efrons_effective_sample_size(
      n = nrow(healthy_data),
      efrons_rms(
        vector2triangle(colMeans(healthy_data), diag_value = 1)
        )
      )
    )

  sick_list <- list(
    data = sick_data,
    jacobian = compute_mu_alpha_jacobian(
      type = 'sick',
      alpha = cov_obj$alpha,
      healthy_dt = healthy_data,
      sick_dt = sick_data,
      d = d,
      linkFun = linkFun),
    expected_value = if(est_mu) {
      triangle2vector(
        linkFun$FUN(
          t = cov_obj$theta,
          a = cov_obj$alpha,
          d = d
        )
      )
    } else colMeans(sick_data),
    solve_Sigma = solve(corrmat_covariance_from_dt(sick_data, est_n = T)),
    df = efrons_effective_sample_size(
      n = nrow(sick_data),
      efrons_rms(
        vector2triangle(colMeans(sick_data), diag_value = 1)
        )
      )
    )
  
  I0 <- compute_gee_raw('I0', healthy_list) + compute_gee_raw('I0', sick_list)
  solve_I0 <- solve(I0)
  I1 <- compute_gee_raw('I1', healthy_list) + compute_gee_raw('I1', sick_list)
  res <- solve_I0 %*% I1 %*% solve_I0
  return(res)
}


compute_fisher_by_grad <- function(cov_obj, sick_dt, linkFun = linkFunctions$multiplicative_identity,
                                   dim_alpha = 1, reg_lambda = 0, reg_p = 2){
  
  sick_dt <- convert_corr_array_to_data_matrix_test(sick_dt)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick_dt))/2
  d <- length(cov_obj$alpha)/p
  
  mu_alpha_jacobian <- 
    jacobian(
      func = function(A) triangle2vector(
        linkFun$FUN(
          t = cov_obj$theta,
          a = A,
          d = d
          )
        ),
      x = cov_obj$alpha
      )
  
  
  g11 <- triangle2vector(
    linkFun$FUN(
      t = cov_obj$theta,
      a = cov_obj$alpha,
      d = length(cov_obj$alpha)/p
    )
  )
  
  residuals <- sick_dt - rep(1, nrow(sick_dt)) %o% g11
  
  Sigma <- triangled_corrmat_covariance(vector2triangle(colMeans(sick_dt), diag_value = 1))/cov_obj$Est_N
  temp_equation <- t(mu_alpha_jacobian) %*% solve(Sigma) %*% t(residuals)
  res <- temp_equation %*% t(temp_equation)
  return(res)
}


compute_fisher_by_hess <- function(cov_obj, sick_dt, linkFun = linkFunctions$multiplicative_identity,
                                   dim_alpha = 1, reg_lambda = 0, reg_p = 2){
  sick_dt <- convert_corr_array_to_data_matrix_test(sick_dt)
  
  Sigma <- triangled_corrmat_covariance(vector2triangle(colMeans(sick_dt), diag_value = 1))/cov_obj$Est_N
  inv_sigma <- solve(Sigma)
  to_deriv <- function(A){
    sum_of_squares(
      alpha = A,
      theta = cov_obj$theta,
      sick_dt = sick_dt,
      inv_sigma = inv_sigma,
      linkFun = linkFun)
  }
  fisher_mat <- hessian(to_deriv, cov_obj$alpha)
  return(fisher_mat)
}


compute_sandwhich_fisher_variance <- function(cov_obj, sick_dt, linkFun = linkFunctions$multiplicative_identity,
                                              dim_alpha = 1, reg_lambda = 0, reg_p = 2){
  grad_fisher <- compute_fisher_by_grad(
    cov_obj = cov_obj,
    sick_dt = sick_dt,
    linkFun = linkFun,
    dim_alpha = dim_alpha,
    reg_lambda = reg_lambda,
    reg_p = reg_p)
  hess_fisher <- compute_fisher_by_hess(
    cov_obj = cov_obj,
    sick_dt = sick_dt,
    linkFun = linkFun,
    dim_alpha = dim_alpha,
    reg_lambda = reg_lambda,
    reg_p = reg_p)
  hess_fisher_solve <- solve(hess_fisher)
  out <- hess_fisher_solve %*% grad_fisher %*% hess_fisher_solve
  return(out)
}


wilks_test <- function(cov_obj, healthy_dt, sick_dt, dim_alpha = 1, linkFun){
  healthy_dt <- convert_corr_array_to_data_matrix_test(healthy_dt)
  sick_dt <- convert_corr_array_to_data_matrix_test(sick_dt)

  #Do a wilks test (chi-square)
  chisq <- -2*( minusloglik(theta = cov_obj$theta,
                            alpha = cov_obj$alpha,
                            linkFun = linkFun,
                            healthy.data = healthy.dat,
                            sick.data = sick.dat,
                            dim_alpha = dim_alpha) -
                  minusloglik(theta = colMeans(rbind(healthy.dat, sick.dat)),
                              alpha = rep(linkFun$NULL_VAL, length(cov_obj$alpha)),
                              linkFun = linkFun,
                              healthy.data = healthy.dat,
                              sick.data = sick.dat,
                              dim_alpha = dim_alpha) )
  
  c("Chisq_val" = chisq, "DF" = length(cov_obj$alpha),
    "Pval" = 1 - pchisq(chisq, length(cov_obj$alpha)))
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

