clean_sick <- function(sick.data, alpha){
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix()
  sick.data/(rep(1, nrow(sick.data)) %*% t(alpha %>% create_alpha_mat() %>% triangle2vector()))
}

minusloglik <- function(theta, alpha, healthy.data, sick.data, effective.N, DET = TRUE){
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
    e20 <- eigen(g20/n.effective_H, symmetric = TRUE)
    U0 <- e20$vectors
    D0 <- e20$values
    
    dist0 <- (healthy.data - rep(1, Nh) %*% t(g10)) %*% U0
  }
  
  g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(alpha)))
  
  g21 <- vector_var_matrix_calc_COR_C(vector2triangle(theta)*create_alpha_mat(alpha))
  if(calc_n) n.effective_D <- compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data), g21)
  e21 <- eigen(g21/n.effective_D, symmetric = TRUE)
  U1 <- e21$vectors
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

cppFunction(paste0(scan("main_work/code/corcalc_c.cpp", what = "character", sep = "\n"), collapse = "\n"))

vector_var_matrix_calc_COR_C <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                         reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i))) - 1
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p)) - 1
  
  pelet <- corcalc_c(MATR, p, m, order_vecti, order_vectj)
  pelet <- pelet + t(pelet) - diag(diag(pelet))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

vector_var_matrix_calc_COR_R <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                         reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i)))
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p))
  
  pelet <- matrix(0, nrow = m, ncol = m)
  for(i1 in 1:m){
    for(j1 in i1:m){
      i <- order_vecti[i1]
      j <- order_vectj[i1]
      k <- order_vecti[j1]
      l <- order_vectj[j1]
      
      MATRij <- MATR[i,j]
      MATRkl <- MATR[k,l]
      MATRik <- MATR[i,k]
      MATRil <- MATR[i,l]
      MATRjk <- MATR[j,k]
      MATRjl <- MATR[j,l]
      
      pelet[i1,j1] <-
        (MATRij*MATRkl/2) * (MATRik^2 + MATRil^2 + MATRjk^2 + MATRjl^2) -
        MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
        MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
        (MATRik*MATRjl + MATRil*MATRjk)
    }
  }
  
  pelet <- pelet + t(pelet) - diag(diag(pelet))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

compute_estimated_N <- function(est, theo){
  x <- diag(theo)
  y <- diag(est)
  return(lm(x ~ 0 + y)$coef)
}

Estimate.Loop <- function(Healthy_List, Sick_List, MaxLoop = 500, Persic = 0.001, method = "BFGS"){
  
  if(class(Healthy_List) == "array") Healthy_List <- cor.matrix_to_norm.matrix(Healthy_List)
  if(class(Sick_List) == "array") Sick_List <- cor.matrix_to_norm.matrix(Sick_List)
  
  healthyMean <- colMeans(Healthy_List)
  sickMean <- colMeans(Sick_List)
  
  N <- c("H" = nrow(Healthy_List), "S" = nrow(Sick_List))
  p <- 0.5 + sqrt(1 + 8*length(sickMean))/2
  
  for.optim <- function(alpha, theta) sum((theta*triangle2vector(create_alpha_mat(alpha)) - sickMean)^2)
  
  Steps <- list()
  temp.theta <- healthyMean
  temp.alpha <- optim(rep(0.8,p), function(A) for.optim(A, temp.theta), method = method)$par
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  
  
  const <- 0
  i <- 1
  distanceA <- 100
  distanceT <- 100
  while((i <= MaxLoop) & (distanceA > Persic) & (distanceT > Persic)){
    
    temp.theta <- rbind(clean_sick(Sick_List, temp.alpha), Healthy_List) %>% colMeans()
    
    temp.alpha <- optim(temp.alpha,function(A) for.optim(A,temp.theta), method = method)$par
    Steps[[i+1]] <- list(theta = temp.theta, alpha = temp.alpha)
    
    distanceA <- vnorm(Steps[[i+1]]$alpha - Steps[[i]]$alpha, sqroot = TRUE)
    distanceT <- vnorm(Steps[[i+1]]$theta - Steps[[i]]$theta, sqroot = TRUE)
    
    if(i%%10==0) cat(paste(i,","))
    i <- i + 1
  }
  cat("\n")
  
  return(list(theta = temp.theta,
              alpha = temp.alpha,
              "Returns" = i,
              Steps = Steps))
}

Estimate.Loop2 <- function(theta0, alpha0, healthy.data, sick.data, T_thresh,
                           max.loop = 50, epsIter = 2*10^(-3), min_reps = 3, method = "Nelder-Mead",
                           epsOptim = 10^(-5), progress = TRUE){
  
  compute_estimated_N_2 <- function(sick.data, theta, alpha, threshold){
    return(min( compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data),
                                    vector_var_matrix_calc_COR_C(vector2triangle(theta)*create_alpha_mat(alpha))),
                threshold ))
  }
  
  temp.theta <- theta0
  temp.alpha <- alpha0
  if( !(is.positive.definite(vector2triangle(theta0)) &
        is.positive.definite(vector2triangle(theta0)*create_alpha_mat(alpha0))) ){
    stop("Initial parameters dont result with positive-definite matrices")
  }
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  healthy_N <- nrow(healthy.data)
  sick_N <- nrow(sick.data)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick.data))/2
  m <- 0.5*p*(p-1)
  
  Steps <- list()
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  log_optim <- list()
  dist <- 100*epsIter
  i <- 1
  
  #Convergence is a matrix wich tells us if the convergence in each iteration is completed
  convergence <- rep(-1, max.loop)
  convergence[1] <- 0
  tnai0 <- FALSE
  
  tt <- Sys.time()
  if(progress) message(paste0("Time of intialization: ", tt, "; Progress: 'Loop, (Time, Convergence, Distance)'"))
  while((i <= max.loop) & (!tnai0)){
    i <- i + 1
    
    effective.N <- compute_estimated_N_2(sick.data, temp.theta, temp.alpha, T_thresh)

    temp.theta <- rbind(healthy.data, clean_sick(sick.data, temp.alpha)) %>% colMeans()
    
    #Optimize Alpha
    optim.alpha <- optim(temp.alpha,
                         function(A) minusloglik(theta = temp.theta,
                                                 alpha = A,
                                                 sick.data = sick.data,
                                                 effective.N = effective.N),
                         method = method,
                         control = list(maxit = min(max(500, i*100), 2000),
                                        reltol = epsOptim))
    convergence[i] <- optim.alpha$convergence
    temp.alpha <- optim.alpha$par
    
    Steps[[i]] <- list(theta = temp.theta, alpha = temp.alpha,
                       convergence = optim.alpha$convergence,
                       Est_N = effective.N)
    log_optim[[i]] <- optim.alpha
    
    #KLAL ATZIRA
    #dist <- max(c(sqrt(vnorm(Steps[[i]]$theta - Steps[[i-1]]$theta)/m),
    #              sqrt(vnorm(Steps[[i]]$alpha - Steps[[i-1]]$alpha)/p)))
    dist <- sqrt(mean((Steps[[i]]$alpha - Steps[[i-1]]$alpha)^2))
    
    tnai0 <- FALSE
    if(i > min_reps) tnai0 <- (dist <= epsIter) & (sum(convergence[i - 0:(min_reps - 1)]) == 0)
    
    if(progress) cat(paste0(i," (",round(as.double.difftime(Sys.time() - tt, units = "secs")), "s, ",
                            convergence[i],", ",round(dist, 5) , "); "))
  }
  if(progress){
    tt <- Sys.time() - tt
    units(tt) <- "secs"
    tt <- as.numeric(tt)
    message(paste0("\nTotal time: ", floor(tt/60), " minutes and ", round(tt %% 60, 1), " seconds."))
  }

  return( list(theta = temp.theta, alpha = temp.alpha, convergence = convergence[1:(min(which(convergence == -1)) - 1)],
               returns = i, Est_N = effective.N, Steps = Steps, Log_Optim = log_optim) )
}
# todo : try updating eigenvector/values every 3 iterations
# todo : check elemntal library, instalation via R
# todo : Build function that does michael and shahar method (4000 comparisons) and compare power
# todo : Add another minus loglik fun
# todo : Add another dimension to alpha