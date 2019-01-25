minusloglik <- function(theta, alpha, healthy.data, sick.data, effective.N, U0, U1, DET = TRUE){
  calcHealth <- !missing(healthy.data)
  
  calc_n <- TRUE
  if(!missing(effective.N)){
    calc_n <- FALSE
    if(calcHealth){
      n.effective_H <- effective.N[1]
      n.effective_D <- effective.N[2]
    } else {
      n.effective_D <- effective.N[1]
    }
  }
  
  Nd <- nrow(sick.data)
  
  if(calcHealth){
    Nh <- nrow(healthy.data)
    
    g10 <- as.matrix(theta)
    
    g20 <- vector_var_matrix_calc_COR_C(vector2triangle(theta))
    if(calc_n) n.effective_H <- compute_estimated_N(cov(healthy.data)*(nrow(healthy.data) - 1)/nrow(healthy.data), g20)
    if(missing(U0)){
      e20 <- eigen(g20/n.effective_H, symmetric = TRUE)
      U0 <- e20$vectors
    } else {
      e20 <- eigen(g20/n.effective_H, symmetric = TRUE, only.values = TRUE)
    }
    D0 <- e20$values
    
    
    dist0 <- (healthy.data - rep(1, Nh) %*% t(g10)) %*% U0
  }
  
  g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(alpha)))
  
  g21 <- vector_var_matrix_calc_COR_C(vector2triangle(theta)*create_alpha_mat(alpha))
  if(calc_n) n.effective_D <- compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data), g21)
  if(missing(U1)){
    e21 <- eigen(g21/n.effective_D, symmetric = TRUE)
    U1 <- e21$vectors
  } else {
    e21 <- eigen(g21/n.effective_D, symmetric = TRUE, only.values = TRUE)
  }
  D1 <- e21$values

  dist1 <- (sick.data - rep(1, Nd) %*% t(g11)) %*% U1
  
  if(calcHealth){
    SSE <- sum(c(diag(dist0 %*% diag(1/D0) %*% t(dist0)), diag(dist1 %*% diag(1/D1) %*% t(dist1))))
    if(DET) SSE <- SSE + Nh*sum(log(D0)) + Nd*sum(log(D1))
  } else {
    SSE <- sum(diag(dist1 %*% diag(1/D1) %*% t(dist1)))
    if(DET) SSE <- SSE + Nd*sum(log(D1))
  }
  
  return(0.5*SSE)
}
