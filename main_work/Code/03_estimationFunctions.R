linkFunctions <- list("Identity" = list(FUN = function(x) x, INV = function(x) x),
                      "Exponent" = list(FUN = function(x) exp(x), INV = function(x) log(x)),
                      "Logarithm" = list(FUN = function(x) log(x), INV = function(x) exp(x)),
                      "Inverse" =  list(FUN = function(x) 1/x, INV = function(x) 1/x),
                      "Benjamini" = list(FUN = function(x) 1/(1 + x), INV = function(x) 1/x - 1) )

clean_sick <- function(sick.data, alpha, dim_alpha = 1){
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  alpha <- matrix(alpha, nc = dim_alpha)
  sick.data/(rep(1, nrow(sick.data)) %*% t(alpha %>% create_alpha_mat(dim_alpha) %>% triangle2vector()))
}

minusloglik <- function(theta, alpha, healthy.data, sick.data, effective.N, dim_alpha = 1, linkFun = function(x) x,
                        U0 = NULL, U1 = NULL, DET = TRUE){
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
    
    g20 <- vector_var_matrix_calc_COR_C(vector2triangle(g10))
    if(calc_n) n.effective_H <- compute_estimated_N(cov(healthy.data)*(nrow(healthy.data) - 1)/nrow(healthy.data), g20)
    if(is.null(U0)){
      e20 <- eigen(g20/n.effective_H, symmetric = TRUE)
      U0 <- e20$vectors
    } else {
      e20 <- eigen(g20/n.effective_H, symmetric = TRUE, only.values = TRUE)
    }
    D0 <- e20$values
    
    
    dist0 <- (healthy.data - rep(1, Nh) %*% t(g10)) %*% U0
  }
  alpha <- matrix(alpha, nrow = dim_alpha)
  g11 <- as.matrix(theta*triangle2vector(create_alpha_mat(linkFun(alpha), dim_alpha = dim_alpha)))
  
  g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11))
  if(calc_n) n.effective_D <- compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data), g21)
  if(is.null(U1)){
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

cppFunction(paste0(scan("main_work/Code/corcalc_c.cpp", what = "character", sep = "\n", quiet = TRUE), collapse = "\n"))

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

Estimate.Loop <- function(iniAlpha, healthy.data, sick.data, linkFun, dim_alpha = 1,
                          MaxLoop = 500, Persic = 0.001, method = "BFGS"){
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  healthyMean <- colMeans(healthy.data)
  sickMean <- colMeans(sick.data)
  
  N <- c("H" = nrow(healthy.data), "S" = nrow(sick.data))
  p <- 0.5 + sqrt(1 + 8*length(sickMean))/2
  if(length(iniAlpha) == 1) iniAlpha <- rep(iniAlpha, p*dim_alpha)
  if(length(iniAlpha) %% p != 0) stop("Initial alpha paramaters aren't of dim [p x dim_alpha]")
  
  for.optim <- function(alpha, theta){
    alpha <- matrix(alpha, nrow = p)
    sum((theta*triangle2vector(create_alpha_mat(linkFun(alpha), dim_alpha = dim_alpha)) - sickMean)^2)
  }
  
  Steps <- list()
  temp.theta <- healthyMean
  temp.alpha <- optim(iniAlpha, function(A) for.optim(A, temp.theta), method = method)$par
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  
  
  const <- 0
  i <- 1
  distanceA <- 100
  distanceT <- 100
  while((i <= MaxLoop) & (distanceA > Persic) & (distanceT > Persic)){
    tt = clean_sick(sick.data, linkFun(temp.alpha), dim_alpha)
    temp.theta <- rbind(clean_sick(sick.data, linkFun(temp.alpha), dim_alpha), healthy.data) %>% colMeans()

    temp.alpha <- optim(temp.alpha, function(A) for.optim(A, temp.theta), method = method)$par
    Steps[[i+1]] <- list(theta = temp.theta, alpha = temp.alpha)
    
    distanceA <- vnorm(Steps[[i+1]]$alpha - Steps[[i]]$alpha, sqroot = TRUE)
    distanceT <- vnorm(Steps[[i+1]]$theta - Steps[[i]]$theta, sqroot = TRUE)
    
    i <- i + 1
  }

  return(list(theta = temp.theta,
              alpha = temp.alpha,
              "Returns" = i,
              Steps = Steps))
}

Estimate.Loop2 <- function(theta0, alpha0, healthy.data, sick.data, T_thresh, linkFun,
                           max.loop = 50, epsIter = 2*10^(-3), min_reps = 3, method = "Nelder-Mead",
                           epsOptim = 10^(-5), updateU = 1, progress = TRUE){
  
  compute_estimated_N_2 <- function(sick.data, theta, alpha, threshold){
    return(min( compute_estimated_N(cov(sick.data)*(nrow(sick.data) - 1)/nrow(sick.data),
                                    vector_var_matrix_calc_COR_C(vector2triangle(theta)*create_alpha_mat(alpha, dim_alpha))),
                threshold ))
  }
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)

  p <- 0.5 + sqrt(1 + 8*ncol(sick.data))/2
  m <- 0.5*p*(p-1)
  dim_alpha <- length(alpha0)/p
  
  if( !(is.positive.definite(vector2triangle(theta0)) &
        is.positive.definite(vector2triangle(theta0)*create_alpha_mat(linkFun(alpha0), dim_alpha))) ){
    stop("Initial parameters dont result with positive-definite matrices")
  }
  
  
  temp.theta <- theta0
  temp.alpha <- alpha0
  
  healthy_N <- nrow(healthy.data)
  sick_N <- nrow(sick.data)
  
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
    
    effective.N <- compute_estimated_N_2(sick.data, temp.theta, linkFun(temp.alpha), T_thresh)
    
    updateU <- 1
    if(updateU != 0) if((i - 2) %% updateU == 0)
      U1 <- eigen(vector_var_matrix_calc_COR_C(vector2triangle(temp.theta) * create_alpha_mat(linkFun(temp.alpha), dim_alpha)),
                  symmetric = T, only.values = F)$vectors
    temp.theta <- rbind(healthy.data, clean_sick(sick.data, linkFun(temp.alpha), dim_alpha)) %>% colMeans()
    
    #Optimize Alpha
    optim.alpha <- optim(temp.alpha,
                         function(A) minusloglik(theta = temp.theta,
                                                 alpha = A, U1 = U1,
                                                 linkFun = linkFun,
                                                 sick.data = sick.data,
                                                 effective.N = effective.N, 
                                                 dim_alpha = dim_alpha),
                         method = method,
                         control = list(maxit = min(max(500, i*100), 2000),
                                        reltol = epsOptim))
    convergence[i] <- optim.alpha$convergence
    temp.alpha <- optim.alpha$par
    
    Steps[[i]] <- list(theta = temp.theta, alpha = temp.alpha,
                       convergence = optim.alpha$convergence,
                       Est_N = effective.N)
    log_optim[[i]] <- optim.alpha
    
    # Stopping rule
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

estimateAlpha <- function(healthy.data, sick.data, T_thresh, linkFun, dim_alpha = 1, updateU = 1, progress = TRUE,
                          INIconfig = list(iniAlpha = 0.8, MaxLoop = 500, Persic = 0.001, method = "BFGS"),
                          FULconfig = list(max.loop = 50, epsIter = 2*10^(-3),
                                           min_reps = 3, method = "Nelder-Mead", epsOptim = 10^(-5))){
  
  if(missing(linkFun)){
    message("No linkFun given in input. Using Identity link function.")
    linkFun <- list(FUN = function(x) x, INV = function(x) x)
  }
  INIconfig <- modifyList(list(iniAlpha = 0.8, MaxLoop = 500, Persic = 0.001, method = "BFGS"), INIconfig)
  FULconfig <- modifyList(list(max.loop = 50, epsIter = 2*10^(-3),
                               min_reps = 3, method = "Nelder-Mead", epsOptim = 10^(-5)),
                          FULconfig)
  

  IID <- Estimate.Loop(healthy.data = healthy.data, sick.data = sick.data, linkFun = linkFun$FUN, dim_alpha = dim_alpha,
                       iniAlpha = linkFun$INV(INIconfig$iniAlpha), MaxLoop = INIconfig$MaxLoop, Persic = INIconfig$Persic, method = INIconfig$method)
  
  COV <- Estimate.Loop2(theta0 = IID$theta, alpha0 = IID$alpha,
                        healthy.data = healthy.data, sick.data = sick.data, T_thresh = T_thresh,
                        linkFun = linkFun$FUN,
                        max.loop = FULconfig$max.loop, epsIter = FULconfig$epsIter, min_reps = FULconfig$min_reps,
                        method = FULconfig$method, epsOptim = FULconfig$epsOptim,
                        updateU = updateU, progress = progress)
  
  return(COV)
}

multipleComparison <- function(healthy.data, sick.data, Tlength,
                               p.adjust.method = p.adjust.methods, test = c("lower", "upper", "both")){
  fisherZ <- function(z) 0.5*log((1 + z)/(1 - z))
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  if(length(p.adjust.method > 1)) p.adjust.method <- p.adjust.method[1]
  if(length(test > 1)) test <- test[3]
  
  H_fisherR <- healthy.data %>% fisherZ() %>% colMeans()
  S_fisherR <- sick.data %>% fisherZ() %>% colMeans()
  
  vars <- (1/nrow(healthy.data) + 1/nrow(sick.data))/(Tlength - 3)
  
  Zvals <- (S_fisherR - H_fisherR)/sqrt(vars)
  if(test == "lower") Pvals <- pnorm(Zvals)
  if(test == "upper") Pvals <- 1 - pnorm(Zvals)
  if(test == "both") Pvals <- 2*pnorm(abs(Zvals), lower.tail = F)
  p.adjust(Pvals, p.adjust.method)
}

# todo : check elemntal library, instalation via R
# todo : Add larger dimensions to alpha