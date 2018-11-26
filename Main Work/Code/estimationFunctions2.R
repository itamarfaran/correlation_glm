clean_sick <- function(sick.data, alpha){
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix()
  sick.data/(rep(1, nrow(sick.data)) %*% t(alpha %>% create_alpha_mat() %>% triangle_to_vector()))
}

minusloglik <- function(theta, alpha, healthy.data = NULL, sick.data, effective.N = NULL, DET = TRUE){
  calcHealth <- length(healthy.data) > 0
  
  calc_n <- TRUE
  if(length(effective.N) > 0){
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
    
    g20 <- vector_var_matrix_calc_COR(vector_to_triangle(theta))
    if(calc_n) n.effective_H <- compute_estimated_N(cov(healthy.data)*(nrow(healthy.data) - 1)/nrow(healthy.data), g20)
    e20 <- eigen(g20/n.effective_H)
    U0 <- e20$vectors
    D0 <- e20$values
    
    dist0 <- (healthy.data - rep(1, Nh) %*% t(g10)) %*% U0
  }
  
  g11 <- as.matrix(theta*triangle_to_vector(create_alpha_mat(alpha)))
  
  g21 <- vector_var_matrix_calc_COR(vector_to_triangle(theta)*create_alpha_mat(alpha))
  if(calc_n) n.effective_D <- compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data), g21)
  e21 <- eigen(g21/n.effective_D)
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
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
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
  
  for.optim <- function(alpha, theta) sum((theta*triangle_to_vector(create_alpha_mat(alpha)) - sickMean)^2)
  
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
  
  return(list(theta = temp.theta,
              alpha = temp.alpha,
              "Returns" = i,
              Steps = Steps))
}

Estimate.Loop2 <- function(theta0, alpha0, healthy.data, sick.data, T_thresh,
                           max.loop = 50, eps = 10^(-3), min_reps = 3, method = "Nelder-Mead",
                           progress = TRUE){
  
  compute_estimated_N_2 <- function(sick.data, theta, alpha, threshold){
    return(min( compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data),
                                    vector_var_matrix_calc_COR(vector_to_triangle(theta)*create_alpha_mat(alpha))),
                threshold ))
  }
  
  temp.theta <- theta0
  temp.alpha <- alpha0
  if( !(is.positive.definite(vector_to_triangle(theta0)) &
        is.positive.definite(vector_to_triangle(theta0)*create_alpha_mat(alpha0))) ){
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

    temp.theta <- rbind(healthy.data, clean_sick(sick.data, temp.alpha)) %>% colMeans()
    
    #Optimize Alpha
    optim.alpha <- optim(temp.alpha,
                         function(A) minusloglik(theta = temp.theta,
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

  return( list(theta = temp.theta, alpha = temp.alpha, convergence = convergence[1:(min(which(convergence == -1)) - 1)],
               returns = i, Est_N = effective.N, Steps = Steps, Log_Optim = log_optim) )
}

loglik_uni <- function(obs, theta, alpha = NULL, Eff.N){
  
  if(length(alpha) == 0){
    meanMat <- vector_to_triangle(theta)
  } else {
    meanMat <- vector_to_triangle(theta)*create_alpha_mat(alpha)
  }
  varMat <- vector_var_matrix_calc_COR(meanMat)/Eff.N
  meanVect <- triangle_to_vector(meanMat)
  
  eigenDec <- eigen(varMat)
  U <- eigenDec$vectors
  D <- eigenDec$values
  
  dists <- t(U) %*% (obs - meanVect)
  
  -0.5*(sum(log(D)) + sum(dists^2/D))
}

loglikgrad_uni <- function(obs, CovObj){
  grad(function(x) loglik_uni(obs = obs, theta = CovObj$theta, alpha = x, Eff.N = CovObj$Est_N), CovObj$alpha) %>%
    (function(x) x %*% t(x))
}

computeBmatr <- function(sickDat, CovObj, ncores = 1){
  Bmatr <- lapply(1:ncol(sickDat), function(j) loglikgrad_uni(sickDat[j,], CovObj)) 
  Bmatr <- Bmatr %>% (function(list){
    L <- length(list)
    pelet <- matrix(0, nrow = nrow(list[[1]]), ncol = ncol(list[[1]]))
    for(i in 1:L) pelet <- pelet + list[[i]]
    pelet
    })
  
  return(Bmatr)
}

ComputeFisher <- function(CovObj, sickDat, method = c("Hess", "Grad")){
  if(class(sickDat) == "array") sickDat <- cor.matrix_to_norm.matrix(sickDat) 
  
  method <- method[1]
  if(method == "Hess") pelet <- hessian(x = CovObj$alpha,
                                        func = function(A) minusloglik(theta = CovObj$theta,
                                                                       alpha = A,
                                                                       sick.data = sickDat,
                                                                       effective.N = CovObj$Est_N))
  if(method == "Grad") pelet <- computeBmatr(sickDat, CovObj)
  
  return(pelet)
}

build_hyp.test <- function(Estimate.Loop2_object, FisherMatr, Real = NULL, test = c("lower", "upper", "two-sided"),
                             qval = 0.05, MH_method = "none", const = 1, effectiveN = NULL){
  
  if(length(test) > 1) test <- test[1]
  obj <- Estimate.Loop2_object
  
  alpha_var_mat <- solve(FisherMatr)
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
  
  A[,5] <- p.adjust(A[,4], method = MH_method)
  A[,6] <- A[,5] < qval
  
  if(length(Real) > 0) A <- cbind(A, Real)
  
  return(list(Results = as.data.frame(A), DF = effectiveN, Test = test, Significance = qval, MH_method = MH_method, SD_origin = alpha_sd, Var_Mat = alpha_var_mat))
}

wilksTest <- function(covObj, healthy.dat, sick.dat){
  if(class(healthy.dat) == "array") healthy.dat <- cor.matrix_to_norm.matrix(healthy.dat)
  if(class(sick.dat) == "array") sick.dat <- cor.matrix_to_norm.matrix(sick.dat)
  
  #Do a wilks test (chi-square)
  chisq <- -2*( minusloglik(theta = covObj$theta,
                            alpha = covObj$alpha,
                            healthy.data = healthy.dat,
                            sick.data = sick.dat) -
                  minusloglik(theta = rbind(healthy.dat, sick.dat) %>% colMeans(),
                              alpha = rep(1, length(covObj$alpha)),
                              healthy.data = healthy.dat,
                              sick.data = sick.dat) )
  
  c("Chisq_val" = chisq, "DF" = length(covObj$alpha),
    "Pval" = 1 - pchisq(chisq, length(covObj$alpha)))
  
}

# minusloglik_onlyalpha <- function(theta, alpha, sick.data, effective.N){
#   Nd <- nrow(sick.data)
#   meanmat <- vector_to_triangle(theta)*create_alpha_mat(alpha)
#   
#   g11 <- as.matrix(triangle_to_vector(meanmat))
#   g21 <- vector_var_matrix_calc_COR(meanmat)/effective.N
#   eigen21 <- eigen(g21)
#   
#   U <- eigen21$vectors
#   D <- eigen21$values
#   dist <- (sick.data - rep(1,nrow(sick.data))%*%t(g11))%*%U
#   
#   return( Nd*sum(log(D)) + sum(diag(dist %*% diag(1/D) %*% t(dist))) )
# }

###

# compute_estimated_N <- function(est, theo){
#   x <- triangle_to_vector(theo, diag = TRUE)
#   y <- triangle_to_vector(est, diag = TRUE)
#   return(lm(x ~ 0 + y)$coef)
# }
