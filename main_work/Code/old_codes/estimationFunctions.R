#Estimate Parameters
Estimate.Loop <- function(Healthy_List, Sick_List, MaxLoop = 500, Persic = 5){
  
  Tnai <- all(abind(Healthy_List, Sick_List, along = 3) %>%
                apply(3, is.positive.definite))
  if(!Tnai) stop("Some matrices not positive definite")
  
  SICK <- Sick_List %>% calculate_mean_matrix() %>% force_symmetry()
  HEALTHY <- Healthy_List %>% calculate_mean_matrix() %>% force_symmetry()
  N <- c(SICK = dim(Sick_List)[3], HEALTHY = dim(Healthy_List)[3])
  p <- dim(Sick_List)[2]
  
  for.optim <- function(ALPHA, THETA) sum(triangle_to_vector(THETA*create_alpha_mat(ALPHA)-SICK)^2)
  
  #Initialize parameters
  Steps <- list()
  temp.theta <- HEALTHY
  temp.alpha <- optim(rep(0.8,p), function(A) for.optim(A,temp.theta), method = "BFGS")$par
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  
  #Begin Loop until Convergence
  const <- 0
  
  i <- 1
  distanceA <- 100
  distanceT <- 100
  while((i<=MaxLoop)&(distanceA>10^(-Persic))&(distanceT>10^(-Persic))){
    temp.theta <- force_symmetry( (N["SICK"] * SICK / create_alpha_mat(temp.alpha) + 
                                     N["HEALTHY"] * HEALTHY) / sum(N) )
    temp.alpha <- optim(temp.alpha,function(A) for.optim(A,temp.theta), method = "BFGS")$par
    Steps[[i+1]] <- list(theta = temp.theta, alpha = temp.alpha)
    
    distanceA <- vnorm(Steps[[i+1]]$alpha-Steps[[i]]$alpha, sqroot = TRUE)
    distanceT <- vnorm(Steps[[i+1]]$theta-Steps[[i]]$theta, sqroot = TRUE)
    
    if(i%%10==0) cat(paste(i,","))
    i <- i + 1
  }
  
  theta.res <- force_positive_definiteness(temp.theta, sensitivity = 0.0001)
  Estimates <- list(theta = theta.res$Matrix,
                    const = theta.res$Alpha,
                    alpha = temp.alpha)
  
  return(list("Estimates" = Estimates, "Returns" = i, Steps = Steps))
}

vector_var_matrix_calc_COR <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  real.cov2 <- function(i, j, k, l) {
    (MATR[i,j]*MATR[k,l]/2) * (MATR[i,k]^2 + MATR[i,l]^2 + MATR[j,k]^2 + MATR[j,l]^2) -
      MATR[i,j]*(MATR[i,k]*MATR[i,l] + MATR[j,k]*MATR[j,l]) -
      MATR[k,l]*(MATR[i,k]*MATR[j,k] + MATR[i,l]*MATR[j,l]) +
      (MATR[i,k]*MATR[j,l] + MATR[i,l]*MATR[j,k])
  }
  
  p <- dim(MATR)[1]
  m <- p*(p-1)/2
  
  v1 <- numeric(0)
  v2 <- numeric(0)
  for(i in 1:(p - 1)){
    v1 <- c(v1, rep(i, p - i))
    v2 <- c(v2, (i + 1):p)
  }
  order_vect <- cbind(v1,v2)
  
  pelet <- matrix(nrow = m, ncol = m)
  for(i in 1:m){
    for(j in i:m){
      indexes <- c(order_vect[i,], order_vect[j,])
      pelet[i,j] <- real.cov2(indexes[1], indexes[2], indexes[3], indexes[4])
      pelet[j,i] <- pelet[i,j]
    }
  }
  
  if((reg_par<0)|(reg_par>1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

###

# compute_estimated_N <- function(est, theo){
#   x <- triangle_to_vector(theo, diag = TRUE)
#   y <- triangle_to_vector(est, diag = TRUE)
#   return(lm(x ~ 0 + y)$coef)
# }

compute_estimated_N <- function(est, theo){
  x <- diag(theo)
  y <- diag(est)
  return(lm(x ~ 0 + y)$coef)
}

minusloglik <- function(theta, alpha, healthy.data, sick.data, effective.N = NULL, DET = TRUE){
  
  calc_n <- TRUE
  if(length(effective.N) == 2){
    calc_n <- FALSE
    n.effective_H <- effective.N[1]
    n.effective_D <- effective.N[2]
  }
  
  Nh <- nrow(healthy.data)
  Nd <- nrow(sick.data)
  
  g20 <- vector_var_matrix_calc_COR(vector_to_triangle(theta))
  if(calc_n) n.effective_H <- compute_estimated_N(cov(healthy.data), g20)
  e20 <- eigen(g20/n.effective_H)
  U0 <- e20$vectors
  D0 <- e20$values
  
  g21 <- vector_var_matrix_calc_COR(vector_to_triangle(theta)*
                                      create_alpha_mat(alpha))
  if(calc_n) n.effective_D <- compute_estimated_N(cov(sick.data), g21)
  e21 <- eigen(g21/n.effective_H)
  U1 <- e21$vectors
  D1 <- e21$values
  
  
  g10 <- as.matrix(theta)
  g11 <- as.matrix(triangle_to_vector(vector_to_triangle(theta)*create_alpha_mat(alpha)))
  
  dist0 <- (healthy.data - rep(1,Nh)%*%t(g10))%*%U0
  dist1 <- (sick.data - rep(1,Nd)%*%t(g11))%*%U1
  
  SSE <- sum(c(diag(dist0%*%diag(1/D0)%*%t(dist0)), diag(dist1%*%diag(1/D1)%*%t(dist1))))
  if(DET) SSE <- SSE + Nh*sum(log(D0)) + Nd*sum(log(D1))
  
  return(SSE)
}

minusloglik_onlyalpha <- function(theta, alpha, sick.data, effective.N){
  Nd <- nrow(sick.data)
  meanmat <- vector_to_triangle(theta)*create_alpha_mat(alpha)
  
  g11 <- as.matrix(triangle_to_vector(meanmat))
  g21 <- vector_var_matrix_calc_COR(meanmat)/effective.N
  eigen21 <- eigen(g21)
  
  U <- eigen21$vectors
  D <- eigen21$values
  dist <- (sick.data - rep(1,nrow(sick.data))%*%t(g11))%*%U
  
  return( Nd*sum(log(D)) + sum(diag(dist %*% diag(1/D) %*% t(dist))) )
}

Estimate.Loop2 <- function(theta0, alpha0, Healthy.ARR, Sick.ARR, T_thresh,
                           max.loop = 50, eps = 10^(-3), min_reps = 3, method = "Nelder-Mead",
                           comp_hess = TRUE, progress = TRUE){
  
  compute_estimated_N_2 <- function(sick.data, theta, alpha, threshold){
    return(min( compute_estimated_N(cov(sick.data), vector_var_matrix_calc_COR(
      vector_to_triangle(theta)*create_alpha_mat(alpha))), threshold ))
  }
  
  temp.theta <- theta0
  temp.alpha <- alpha0
  if( !(is.positive.definite(vector_to_triangle(theta0)) &
        is.positive.definite(vector_to_triangle(theta0)*create_alpha_mat(alpha0))) ){
    stop("Initial parameters dont result with positive-definite matrices")
  }
  
  healthy.data <- cor.matrix_to_norm.matrix(Healthy.ARR)
  sick.data <- cor.matrix_to_norm.matrix(Sick.ARR)
  healthy_N <- nrow(healthy.data)
  sick_N <- nrow(sick.data)
  p <- dim(Healthy.ARR)[1]
  m <- 0.5*p*(p-1)
  
  Steps <- list()
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  log_optim <- list()
  dist <- 100*eps
  i <- 1
  
  #Convergence is a matrix wich tells us if the convergence in each iteration is completed
  convergence <- rep(-1, max.loop)
  convergence[1] <- 0
  tnai0 <- FALSE
  
  if(progress) cat("\n Progress: 'Loop, (Convergence, Distance);'\n")
  while((i <= max.loop) & (!tnai0)){
    i <- i + 1
    
    effective.N <- compute_estimated_N_2(sick.data, temp.theta, temp.alpha, T_thresh)
    #effective.N <- Tlength
    
    clean_sick <- sick.data/(rep(1, sick_N)
                             %*% t(temp.alpha %>% create_alpha_mat() %>% triangle_to_vector()))
    
    temp.theta <- rbind(healthy.data, clean_sick) %>% colMeans()
    
    #Optimize Alpha
    optim.alpha <- optim(temp.alpha,
                         function(A) minusloglik_onlyalpha(theta = temp.theta,
                                                           alpha = A,
                                                           sick.data = sick.data,
                                                           effective.N = effective.N),
                         method = method,
                         control = list(maxit = min(max(500, i*100), 2000)))
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
    if(i > min_reps) tnai0 <- (dist <= eps) & (sum(convergence[i - 0:(min_reps - 1)]) == 0)
    
    if(progress) cat(paste0(i," (",convergence[i],", ",round(dist, 5) , "); "))
  }
  if(comp_hess){
    cat("\n \n 'while' loop finished; Computing Hessian\n")
    
    hess <- hessian(x = temp.alpha,
                    func = function(A) minusloglik_onlyalpha(theta = temp.theta,
                                                             alpha = A,
                                                             sick.data = sick.data,
                                                             effective.N = effective.N))
  } else {hess <- NULL}
  
  return( list(theta = temp.theta, alpha = temp.alpha, Hess = hess, convergence = convergence[1:(min(which(convergence == -1)) - 1)],
               returns = i, Est_N = effective.N, Steps = Steps, Log_Optim = log_optim) )
}

# build_hyp.test <- function(Estimate.Loop2_object, Real = NULL, test = c("lower", "upper", "two-sided"),
#                            qval = 0.05, method = "none", const = 1, effectiveN = NULL){
build_hyp.test <- function(Estimate.Loop2_object, Real = NULL, test = c("lower", "upper", "two-sided"),
                             qval = 0.05, method = "none", const = 1, effectiveN = NULL){
    if(length(test) > 1) test <- test[1]
  obj <- Estimate.Loop2_object
  
  # alpha_var_mat <- solve(obj$Hess)
  alpha_var_mat <- solve(obj$Hess)
  alpha_sd <- sqrt(diag(alpha_var_mat))
  
  A <- matrix(nrow = length(obj$alpha), ncol = 6)
  colnames(A) <- c("Est.", "Std.", "Z-val", "P-val", "Adj.P-val", "Reject_H0")
  
  dist_fun <- function(q) pnorm(q)
  if(length(effectiveN) > 0){
    colnames(A)[3] <- "T-val"
    dist_fun <- function(q) pt(q, effectiveN)
  }
  
  A[,1] <- obj$alpha
  A[,2] <- const*alpha_sd
  A[,3] <- (A[,1] - 1)/A[,2]
  
  if(test == "lower"){ A[,4] <- round(dist_fun(A[,3]),5)
  }else if(test == "upper"){ A[,4] <- round(1 - dist_fun(A[,3]),5)
  }else{A[,4] <- round(2*(1-dist_fun(abs(A[,3]))),5)}
  
  A[,5] <- p.adjust(A[,4], method = method)
  A[,6] <- A[,5] < qval
  
  if(length(Real) > 0) A <- cbind(A, Real)
  
  return(list(Results = as.data.frame(A), DF = effectiveN, Test = test, Significance = qval, method = method, SD_origin = alpha_sd, Var_Mat = alpha_var_mat))
}

build_correction_mat <- function(Estimate.Loop2_object, Sick.ARR){
  sick.data <- cor.matrix_to_norm.matrix(Sick.ARR)
  theta <- Estimate.Loop2_object$theta
  alpha <- Estimate.Loop2_object$alpha
  k <- nrow(sick.data)
  
  forDeriv1 <- function(X){
    g21 <- vector_var_matrix_calc_COR(vector_to_triangle(theta) * create_alpha_mat(X))
    sum(log(eigen(g21, only.values = TRUE)$values))
  }
  psi <- grad(forDeriv1, alpha)
  
  forDeriv2 <- function(X, i){
    x <- sick.data[i,]
    meanmat <- vector_to_triangle(theta)*create_alpha_mat(X)
    t(x - triangle_to_vector(meanmat)) %*% solve(vector_var_matrix_calc_COR(meanmat)) %*% (x - triangle_to_vector(meanmat))
  } 
  
  xi <- t(sapply(1:k, function(i) grad(forDeriv2, alpha, i = i)))
  
  p <- length(psi)
  
  sumMat <- k*((psi) %*% t(psi))
  for(i in 1:k){
    psixi <- psi %*% t(xi[i,])
    sumMat <- sumMat + psixi + t(psixi) + (xi[i,] %*% t(xi[i,]))
  }
  
  sumMat <- sumMat/(4*k)
  HessSolve <- -solve(Estimate.Loop2_object$Hess)
  return(list(Correction = sumMat,
              Corrected = HessSolve %*% sumMat %*% HessSolve))
}

Estimate.Loop2_old <- function(theta0, alpha0, healthy.data, sick.data, T_thresh,
                               linkFun = linkFunctions$multiplicative_identity,
                               var_weights = c(1, 0, 0), sigma = 1, reg_lambda = 0, reg_p = 2,
                               nonpositive = "Stop", max.loop = 50, epsIter = 2*10^(-3), min_reps = 3,
                               method = "Nelder-Mead", epsOptim = 10^(-5), updateU = 1, progress = TRUE){
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick.data))/2
  m <- 0.5*p*(p-1)
  dim_alpha <- length(alpha0)/p
  if(dim_alpha %% 1 != 0) stop("alpha0 not multiplicative of p")
  
  if( !(is.positive.definite(vector2triangle(theta0)) &
        is.positive.definite(linkFun$FUN(t = theta0, a = alpha0, d = dim_alpha))) ){
    stop("Initial parameters dont result with positive-definite matrices")
  }
  
  temp.theta <- theta0
  temp.alpha <- alpha0
  
  healthy_N <- nrow(healthy.data)
  sick_N <- nrow(sick.data)
  
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
      vector_var_matrix_calc_COR_C(g11, nonpositive = "Ignore"))
    effective.N <- min(effective.N, T_thresh)
    
    temp.theta <- colMeans(rbind(
      linkFun$CLEAN(dt = sick.data, a = temp.alpha, d = dim_alpha),
      healthy.data
    ))
    
    if(updateU == 0){
      U1 <- NULL
    } else if((i - 2) %% updateU == 0){
      U1 <- eigen(vector_var_matrix_calc_COR_C(g11, nonpositive = nonpositive),
                  symmetric = T, only.values = F)$vectors
    }
    
    #Optimize Alpha
    optim.alpha <- optim(
      temp.alpha,
      function(A) {
        minusloglik(theta = temp.theta,
                    alpha = A, U1 = U1,
                    linkFun = linkFun,
                    sick.data = sick.data,
                    effective.N = effective.N,
                    dim_alpha = dim_alpha,
                    nonpositive = nonpositive,
                    var_weights = var_weights,
                    sigma = sigma) +
          reg_lambda*sum((A - linkFun$NULL_VAL)^reg_p)
      },
      method = method,
      control = list(maxit = min(max(500, i*100), 2000),
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
  
  return( list(theta = temp.theta, alpha = temp.alpha,
               convergence = convergence[1:(min(which(convergence == -1)) - 1)],
               returns = i, Est_N = effective.N,
               Steps = Steps, Log_Optim = log_optim) )
}
