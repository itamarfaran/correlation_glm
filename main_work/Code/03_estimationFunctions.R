compute_estimated_N <- function(est, theo){
  x <- diag(theo)
  y <- diag(est)
  return(lm(x ~ 0 + y)$coef)
}


cppFunction(paste0(scan("main_work/Code/corcalc_c.cpp", what = "character", sep = "\n", quiet = TRUE), collapse = "\n"))


vector_var_matrix_calc_COR_C <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                         reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {
      MATR <- as.matrix(Matrix::nearPD(MATR, corr = TRUE, doSym = TRUE)$mat)
    }
    else if(nonpositive != "Ignore") {
      stop("MATR not positive definite")
    }
  }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i))) - 1
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p)) - 1
  
  output <- corcalc_c(MATR, p, m, order_vecti, order_vectj)
  output <- output + t(output) - diag(diag(output))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) output <- (1 - reg_par)*output + reg_par*diag(diag(output))
  
  return(output)
}

vector_var_matrix_calc_COR_R <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                         reg_par = 0){
  
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {
      MATR <- as.matrix(Matrix::nearPD(MATR, corr = TRUE, doSym = TRUE)$mat)
    }
    else if(nonpositive != "Ignore") {
      stop("MATR not positive definite")
    }
  }
  
  p <- nrow(MATR)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i)))
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p))
  
  output <- matrix(0, nrow = m, ncol = m)
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
      
      output[i1,j1] <-
        (MATRij*MATRkl/2) * (MATRik^2 + MATRil^2 + MATRjk^2 + MATRjl^2) -
        MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
        MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
        (MATRik*MATRjl + MATRil*MATRjk)
    }
  }
  
  output <- output + t(output) - diag(diag(output))
  
  if((reg_par < 0) | (reg_par > 1)) warning("Regularization Parameter not between 0,1")
  if(reg_par != 0) output <- (1 - reg_par)*output + reg_par*diag(diag(output))
  
  return(output)
}


minusloglik <- function(theta, alpha, healthy.data, sick.data, effective.N,
                        linkFun = linkFunctions$multiplicative_identity,
                        var_weights = c(1, 0, 0), sigma = 1, dim_alpha = 1,
                        nonpositive = "Stop", U0 = NULL, U1 = NULL, DET = TRUE){
  calcHealth <- !missing(healthy.data)
  length(var_weights) <- 3; var_weights[is.na(var_weights)] <- 0
  var_weights <- var_weights/sum(var_weights)
  
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
  if(calcHealth){
    Nh <- nrow(healthy.data)
    g10 <- as.matrix(theta)
    
    g20 <- vector_var_matrix_calc_COR_C(vector2triangle(g10), nonpositive = nonpositive)
    g20_adjusted <- matrix(0, nrow = ncol(healthy.data), ncol = ncol(healthy.data))
    if(var_weights[1] > 0) g20_adjusted <- g20_adjusted + var_weights[1]*g20
    if(var_weights[2] > 0) g20_adjusted <- g20_adjusted +
      var_weights[2]*vector_var_matrix_calc_COR_C(vector2triangle(colMeans(healthy.data)))
    if(var_weights[3] > 0) g20_adjusted <- g20_adjusted + var_weights[3]*sigma*diag(ncol(healthy.data))
    
    if(calc_n) n.effective_H <- compute_estimated_N(cov(healthy.data)*(Nh - 1)/Nh, g20)
    if(is.null(U0)){
      e20 <- eigen(g20/n.effective_H, symmetric = TRUE, only.values = FALSE)
      U0 <- e20$vectors
    } else {
      e20 <- eigen(g20/n.effective_H, symmetric = TRUE, only.values = TRUE)
    }
    D0 <- e20$values
    
    dist0 <- (healthy.data - rep(1, Nh) %*% t(g10)) %*% U0
  }
  
  Nd <- nrow(sick.data)
  g11 <- as.matrix(triangle2vector(linkFun$FUN(t = theta, a = alpha, d = dim_alpha)))

  g21 <- vector_var_matrix_calc_COR_C(vector2triangle(g11), nonpositive = nonpositive)
  g21_adjusted <- matrix(0, nrow = ncol(sick.data), ncol = ncol(sick.data))
  if(var_weights[1] > 0) g21_adjusted <- g21_adjusted + var_weights[1]*g21
  if(var_weights[2] > 0) g21_adjusted <- g21_adjusted +
    var_weights[2]*vector_var_matrix_calc_COR_C(vector2triangle(colMeans(sick.data)))
  if(var_weights[3] > 0) g21_adjusted <- g21_adjusted + var_weights[3]*sigma*diag(ncol(sick.data))
  
  if(calc_n) n.effective_D <- compute_estimated_N(cov(sick.data)*(Nd - 1)/Nd, g21)
  if(is.null(U1)){
    e21 <- eigen(g21_adjusted/n.effective_D, symmetric = TRUE, only.values = FALSE)
    U1 <- e21$vectors
  } else {
    e21 <- eigen(g21_adjusted/n.effective_D, symmetric = TRUE, only.values = TRUE)
  }
  D1 <- e21$values
  
  dist1 <- (sick.data - rep(1, Nd) %*% t(g11)) %*% U1
  
  SSE <- 0
  if(calcHealth){
    SSE <- SSE + sum(diag(dist0 %*% diag(1/D0) %*% t(dist0)))
    if(DET) SSE <- SSE + Nh*sum(log(D0))
  }
  SSE <- SSE + sum(diag(dist1 %*% diag(1/D1) %*% t(dist1)))
  if(DET) SSE <- SSE + Nd*sum(log(D1))
    
  return(0.5*SSE)
}

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
  if(length(iniAlpha) == 1) iniAlpha <- rep(iniAlpha, p*dim_alpha)
  if(length(iniAlpha) %% p != 0) stop("Initial alpha paramaters aren't of dim [p x dim_alpha]")
  
  for.optim <- function(alpha, theta){
    alpha <- matrix(alpha, nrow = p)
    sum(
      (triangle2vector(linkFun$FUN(t = theta, a = alpha, d = dim_alpha)) - sickMean)^2
      ) +
      reg_lambda*sum((alpha - linkFun$NULL_VAL)^reg_p)
    }
  
  Steps <- list()
  temp.theta <- healthyMean
  temp.alpha <- optim(iniAlpha, function(A) for.optim(A, temp.theta), method = method)$par
  Steps[[1]] <- list(theta = temp.theta, alpha = temp.alpha)
  
  i <- 1
  distanceA <- 100
  distanceT <- 100
  for(i in 1:MaxLoop){
    temp.theta <- colMeans(rbind(
      linkFun$CLEAN(dt = sick.data, a = temp.alpha, d = dim_alpha),
      healthy.data
      ))

    temp.alpha <- optim(temp.alpha, function(A) for.optim(A, temp.theta), method = method)$par
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

estimateAlpha <- function(healthy.data, sick.data, T_thresh = Inf,
                          linkFun = linkFunctions$multiplicative_identity,
                          dim_alpha = 1, var_weights = c(1, 0, 0), sigma = 1, reg_lambda = 0, reg_p = 2,
                          updateU = 1, progress = TRUE,
                          INIconfig = list(iniAlpha = 0.8, MaxLoop = 500, Persic = 0.001, method = "BFGS"),
                          FULconfig = list(max.loop = 50, epsIter = 2*10^(-3),
                                           min_reps = 3, nonpositive = "Ignore",
                                           method = "Nelder-Mead", epsOptim = 10^(-5))){
  
  INIconfig <- modifyList(list(iniAlpha = 0.8, MaxLoop = 500, Persic = 0.001, method = "BFGS"), INIconfig)
  FULconfig <- modifyList(list(max.loop = 50, epsIter = 2*10^(-3), min_reps = 3, nonpositive = "Ignore",
                               method = "Nelder-Mead", epsOptim = 10^(-5)), FULconfig)
  
  IID <- Estimate.Loop(iniAlpha = linkFun$INV(INIconfig$iniAlpha),
                       healthy.data = healthy.data, sick.data = sick.data,
                       linkFun = linkFun, reg_lambda = reg_lambda, reg_p = reg_p, dim_alpha = dim_alpha,
                       MaxLoop = INIconfig$MaxLoop, Persic = INIconfig$Persic, method = INIconfig$method)
  
  COV <- Estimate.Loop2(theta0 = IID$theta, alpha0 = IID$alpha,
                        healthy.data = healthy.data, sick.data = sick.data, T_thresh = T_thresh,
                        linkFun = linkFun, var_weights = var_weights,
                        sigma = sigma, reg_lambda = reg_lambda, reg_p = reg_p,
                        nonpositive = FULconfig$nonpositive, max.loop = FULconfig$max.loop,
                        epsIter = FULconfig$epsIter, min_reps = FULconfig$min_reps,
                        method = FULconfig$method, epsOptim = FULconfig$epsOptim,
                        updateU = updateU, progress = progress)
  
  
  return(COV)
}

# todo : check elemntal library, instalation via R
# todo : Add larger dimensions to alpha