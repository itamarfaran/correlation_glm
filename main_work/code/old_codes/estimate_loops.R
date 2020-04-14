Estimate.Loop <- function(iniAlpha, healthy.data, sick.data,
                          linkFun = linkFunctions$multiplicative_identity,
                          reg_lambda = 0, reg_p = 2, dim_alpha = 1,
                          MaxLoop = 500, Persic = 0.001, method = "BFGS"){
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  healthyMean <- colMeans(healthy.data)
  sickMean <- colMeans(sick.data)
  
  N <- c("H" = nrow(healthy.data), "S" = nrow(sick.data))
  p <- 0.5 + sqrt(1 + 8*length(sickMean))/2
  m <- 0.5*p*(p-1)
  
  if(length(iniAlpha) == 1) iniAlpha <- rep(iniAlpha, p*dim_alpha)
  if(length(iniAlpha) %% p != 0) stop("Initial alpha paramaters aren't of dim [p x dim_alpha]")
  
  Steps <- list()
  temp.theta <- healthyMean
  temp.alpha <- optim(par = iniAlpha, fn = sum_of_squares, 
                      theta = temp.theta, sick.data = sick.data, inv_sigma = diag(m),
                      linkFun = linkFun, dim_alpha = dim_alpha,
                      reg_lambda = reg_lambda, reg_p = reg_p,
                      method = method)$par
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  
  i <- 1
  distanceA <- 100
  distanceT <- 100
  for(i in 1:MaxLoop){
    
    temp.theta <- theta_of_alpha(
      alpha = temp.alpha,
      healthy_dt = healthy.data,
      sick_dt = sick.data,
      linkFun = linkFun,
      d = dim_alpha) 
    
    temp.alpha <- optim(par = temp.alpha, fn = sum_of_squares, 
                        theta = temp.theta, sick.data = sick.data, inv_sigma = diag(m),
                        linkFun = linkFun, dim_alpha = dim_alpha,
                        reg_lambda = reg_lambda, reg_p = reg_p,
                        method = method)$par
    Steps[[i+1]] <- list(theta = temp.theta, alpha = temp.alpha)
    
    distanceA <- vnorm(Steps[[i+1]]$alpha - Steps[[i]]$alpha, sqrt = TRUE)
    distanceT <- vnorm(Steps[[i+1]]$theta - Steps[[i]]$theta, sqrt = TRUE)
    if((distanceA < Persic) & (distanceT < Persic)) break()
  }
  
  return(list(theta = temp.theta,
              alpha = temp.alpha,
              "Returns" = i,
              Steps = Steps))
}

Estimate.Loop2 <- function(theta0, alpha0, healthy.data, sick.data, T_thresh,
                           linkFun = linkFunctions$multiplicative_identity,
                           var_weights = c(1, 0), sigma = 1, reg_lambda = 0, reg_p = 2,
                           nonpositive = "Stop", max.loop = 50, epsIter = 2*10^(-4), min_reps = 3,
                           method = "BFGS", epsOptim = 10^(-5), updateU = 1, progress = TRUE){
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick.data))/2
  m <- 0.5*p*(p-1)
  dim_alpha <- length(alpha0)/p
  if(dim_alpha %% 1 != 0) stop("alpha0 not multiplicative of p")
  
  if(!(
    is.positive.definite(vector2triangle(theta0)) &
    is.positive.definite(linkFun$FUN(t = theta0, a = alpha0, d = dim_alpha))
  ) ){
    stop("Initial parameters dont result with positive-definite matrices")
  }
  
  temp.theta <- theta0
  temp.alpha <- alpha0
  
  healthy_N <- nrow(healthy.data)
  sick_N <- nrow(sick.data)
  
  g12 <- vector_var_matrix_calc_COR_C(
    vector2triangle(
      colMeans(
        sick.data
      ), diag = FALSE
    )
  )
  
  g12_reg <- var_weights[1]*g12 + var_weights[2]*diag(m)*sigma
  solve_g12_reg <- solve(g12_reg)
  
  Steps <- list()
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  log_optim <- list()
  dist <- 1000*epsIter
  
  #Convergence is a matrix wich tells us if the convergence in each iteration is completed
  convergence <- rep(-1, max.loop)
  convergence[1] <- 0
  tnai0 <- FALSE
  
  tt <- Sys.time()
  if(progress) message(paste0("Time of intialization: ", tt, "; Progress: 'Loop, (Time, Convergence, Distance)'"))
  for(i in 2:max.loop){
    g11 <- linkFun$FUN(t = temp.theta, a = temp.alpha, d = dim_alpha)
    
    effective.N <- compute_estimated_N(
      ((nrow(sick.data) - 1)/nrow(sick.data)) * cov(sick.data),
      g12)
    effective.N <- min(effective.N, T_thresh)
    
    # residuals <- sick.data - (rep(1, sick_N) %o% triangle2vector(g11))
    # theo_sigma <- t(residuals) %*% residuals
    # 
    # effective.N <- compute_estimated_N(
    #   ((nrow(sick.data) - 1)/nrow(sick.data)) * cov(sick.data),
    #   theo_sigma)
    # effective.N <- min(effective.N, T_thresh)
    
    temp.theta <- theta_of_alpha(
      alpha = temp.alpha,
      healthy_dt = healthy.data,
      sick_dt = sick.data,
      linkFun = linkFun,
      d = dim_alpha) 
    
    optim.alpha <- optim(
      par = temp.alpha,
      fn = sum_of_squares,
      theta = temp.theta, sick.data = sick.data, inv_sigma = effective.N*solve_g12_reg,
      linkFun = linkFun, dim_alpha = dim_alpha,
      reg_lambda = reg_lambda, reg_p = reg_p,
      method = method, control = list(
        maxit = min(max(500, i*100), 2000),
        reltol = epsOptim)
    )
    
    convergence[i] <- optim.alpha$convergence
    temp.alpha <- optim.alpha$par
    
    Steps[[i]] <- list(theta = temp.theta, alpha = temp.alpha,
                       convergence = optim.alpha$convergence,
                       Est_N = effective.N)
    log_optim[[i]] <- optim.alpha
    
    
    if(progress) cat(paste0(i," (",round(as.double.difftime(Sys.time() - tt, units = "secs")), "s, ",
                            convergence[i],", ",round(dist, 5) , "); "))
    
    # Stopping rule
    dist <- sqrt(mean((Steps[[i]]$alpha - Steps[[i-1]]$alpha)^2))
    tnai0 <- FALSE
    if(i > min_reps) tnai0 <- (dist <= epsIter) & (sum(convergence[i - 0:(min_reps - 1)]) == 0)
    if(tnai0) break()
  }
  if(progress){
    tt <- Sys.time() - tt
    units(tt) <- "secs"
    tt <- as.numeric(tt)
    message(paste0("\nTotal time: ", floor(tt/60), " minutes and ", round(tt %% 60, 1), " seconds."))
  }
  
  suppressWarnings({
    max_convergence <- min(
      (min(which(convergence == -1)) - 1),
      max.loop
    )
  })
  
  return( list(theta = temp.theta, alpha = temp.alpha, Est_N = effective.N, linkFun = linkFun,
               convergence = convergence[1:max_convergence],
               returns = i, Steps = Steps, Log_Optim = log_optim) )
}

