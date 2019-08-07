### todo: need methodology for simulating alphas

linkFunctions <- list(
  "multiplicative_identity" = list(
    FUN = function(t, a, d) vector2triangle(t)*create_alpha_mat(a, d),
    INV = function(a) a,
    CLEAN = function(dt, a, d)
      return(dt / (rep(1, nrow(dt)) %*% t(a %>% create_alpha_mat(d) %>% triangle2vector()))),
    NULL_VAL = 1
  ),
  "additive_identity" = list(
    FUN = function(t, a, d) vector2triangle(t)*(1 + create_alpha_mat(a, d)),
    INV = function(a) a,
    CLEAN = function(dt, a, d)
      return(dt / (1 + rep(1, nrow(dt)) %*% t(a %>% create_alpha_mat(d) %>% triangle2vector()))),
    NULL_VAL = 0
  ),
  "additive_quotent" = list(
    FUN = function(t, a, d) vector2triangle(t)/(1 + create_additive_alpha_mat(a)),
    INV = function(a) a,
    CLEAN = function(dt, a, d)
      return(dt * (1 + rep(1, nrow(dt)) %*% t(a %>% create_additive_alpha_mat() %>% triangle2vector()))),
    NULL_VAL = 0
  )
)


build_parameters <- function(p, percent_alpha, range_alpha, dim_alpha = 1, loc_scale = c(0,1), seed){
  #Build Real Sigma and Theta
  if(!missing(seed)) set.seed(seed[1])
  temp <- matrix(rnorm(2*p^2, loc_scale[1], loc_scale[2]), nrow = 2*p)
  temp <- t(temp)%*%temp
  real.theta <- force_symmetry(cov2cor(temp))
  
  #Generate Variances and create Variance matrix parameters
  if(!missing(seed)) set.seed(seed[2])
  varss <- rlogis(p)^2
  real.sigma <- sqrt(diag(varss)) %*% real.theta %*% sqrt(diag(varss))
  
  #Build Real Alpha
  sum_alpha <- rep(1, p)
  sum_alpha[sample(p, floor(percent_alpha * p))] <- runif(floor(percent_alpha * p), range_alpha[1], range_alpha[2])
  alpha <- matrix(runif(p*dim_alpha), nr = p)
  alpha <- apply(alpha, 1, function(x) x/sum(x))
  if(dim_alpha > 1) alpha <- t(alpha)
  alpha <- alpha * (sum_alpha %o% rep(1, dim_alpha))
  
  return(list(Corr.mat = real.theta, Cov.mat = real.sigma, Alpha = alpha))
}

create_correlation_matrices <- function(real_corr, real_var, sample_size, df = 0, AR = NULL, MA = NULL, 
                                        seed.control, ncores = 1, silent = FALSE){
  if(missing(real_var)){
    if(!is.positive.definite(real_corr)) stop("real_corr not positive definite")
    p <- nrow(real_corr)
    var_scale <- runif(p,10,100)
    Dhalf <- sqrt(diag(var_scale))
    real_var <- force_symmetry(Dhalf%*%real_corr%*%Dhalf)
  } else {
    if(!missing(real_corr)) warning("Both real_corr and real_var provided; Using real_var")
  }
  if(!is.positive.definite(real_var)) stop("real_var not positive definite")
  if(is.null(AR) & is.null(MA)){
    if(!missing(seed.control)){
      if(df >= p) {
        set.seed(seed.control)
      } else {
        stop("seed.control not available for rWishart2")
      }
    }
    pelet_matrices <- rWishart2(sample_size, df, real_var, ncores = ncores)
  } else {
    if(!missing(seed.control)) stop("seed.control not available for rWishart_ARMA")
    pelet_matrices <- rWishart_ARMA(sample_size, df, real_var, AR = AR, MA = MA, silent = silent, ncores = ncores)
  }

  return( simplify2array(mclapply(1:sample_size, function(b) force_symmetry(cov2cor(pelet_matrices[,,b])),
                                  mc.cores = ncores )) )
}

rWishart2 <- function(n = 1, df, Sigma, ncores = 1){
  p <- ncol(Sigma)
  
  if(df >= p) return(rWishart(n, df, Sigma))
  warning("Wishart degrees of freedom lower than matrix dimension.")
  
  rawFun <- function(k){
    matrices <- rmvnorm(n = df, sigma = Sigma)
    return(t(matrices) %*% matrices)
  }
  
  simplify2array(mclapply(1:n, rawFun, mc.cores = ncores))
}

rWishart_ARMA <- function(n = 1, df, Sigma, AR = NULL, MA = NULL, silent = FALSE, ncores = 1){
  p <- ncol(Sigma)
  
  if(is.null(MA) & is.null(AR)) return(rWishart2(n = n, df = df, Sigma = Sigma))
  if(df < p) warning("Wishart degrees of freedom lower than matrix dimension.")
  
  maxAR <- maxMA <- NULL
  
  
  if(!is.null(AR)){
    if(!checkInv(AR)) stop("AR process not stationary.")
    maxAR <- length(AR)
  }
  if(!is.null(MA)){
    if(!checkInv(MA)) stop("MA process not invertable.")
    maxMA <- length(MA)
  }

  rawFun <- function(k){
    NormMatrix <- rmvnorm(n = df, sigma = Sigma)
    NormMatrix_ARMA <- NormMatrix
    
    for(i in 2:df){
      if(!is.null(AR)){
        arlag <- min(maxAR, i - 1)
        NormMatrix_ARMA[i,] <- NormMatrix_ARMA[i,] + sumvector(NormMatrix_ARMA, (i - arlag):(i - 1), rev(AR[1:arlag]))
        # X[i] = epsilon[i] + X[i-1] + ...
      }
      if(!is.null(MA)){
        malag <- min(maxMA, i - 1)
        NormMatrix_ARMA[i,] <- NormMatrix_ARMA[i,] + sumvector(NormMatrix, (i - malag):(i - 1) , rev(MA[1:malag]))
        # X[i] = epsilon[i] + X[i-1] + ... + epsilon[i-1] + ...
      }
    }
    return(t(NormMatrix_ARMA) %*% NormMatrix_ARMA)
  }
  
  return(simplify2array(mclapply(1:n, rawFun, mc.cores = ncores)))
}

createSamples <- function(B = 1, nH, nS, p, Tlength, percent_alpha, range_alpha, dim_alpha = 1, loc_scale = c(0,1),
                          linkFun = linkFunctions$multiplicative_identity,
                          ARsick = NULL, ARhealth = NULL, MAsick = NULL, MAhealth = NULL, seed = NULL, ncores = 1){
  
  parameters <- build_parameters(p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, dim_alpha = dim_alpha,
                                 loc_scale = loc_scale, seed = seed)
  real.theta <- parameters$Corr.mat
  real.sigma <- parameters$Cov.mat
  alpha <- parameters$Alpha
  g21 <- linkFun$FUN(t = triangle2vector(real.theta), a = alpha, d = dim_alpha)
  
  if(B == 1) return(list(healthy = create_correlation_matrices(real_corr = real.theta, sample_size = nH,
                                                               df = Tlength, AR = ARhealth, MA = MAhealth, ncores = ncores),
                           sick = create_correlation_matrices(real_corr = g21, sample_size = nS,
                                                              df = Tlength, AR = ARsick, MA = MAsick, ncores = ncores),
                         real.theta = real.theta, real.sigma = real.sigma, alpha = alpha))
  
  pelet <- list(real.theta = real.theta, real.sigma = real.sigma, alpha = alpha, samples = list())
  
  rawFun <- function(b){
    list(healthy = create_correlation_matrices(real_corr = real.theta, sample_size = nH, df = Tlength,
                                               AR = ARsick, MA = MAsick, silent = TRUE, ncores = 1),
         sick = create_correlation_matrices(real_corr = g21, sample_size = nS, df = Tlength,
                                            AR = ARsick, MA = MAsick, silent = TRUE, ncores = 1))
  }
  
  pelet$samples <- mclapply(1:B, rawFun, mc.cores = ncores)
  
  return(pelet)
}

rWishart_ARMA2 <- function(n = 1, df, Sigma, AR = NULL, MA = NULL, ncores = 1){
  p <- ncol(Sigma)
  
  if(is.null(MA) & is.null(AR)) return(rWishart2(n = n, df = df, Sigma = Sigma))
  if(df < p) warning("Wishart degrees of freedom lower than matrix dimension.")
  
  maxAR <- maxMA <- NULL
  
  if(!is.null(AR)){
    if(!checkInv(AR)) stop("AR process not stationary.")
    maxAR <- length(AR)
    
  } 
  if(!is.null(MA)){
    if(!checkInv(MA)) stop("MA process not invertable.")
    maxMA <- length(MA)
  }
  
  rawFun <- function(k){
    
    matrices <- rmvnorm(n = df, sigma = Sigma)
    matrices <- sapply(1:df, function(i) matrices[i,] %*% t(matrices[i,]), simplify = "array")
    
    centralMatrices <- array(0, c(p, p, df))
    corelMatrics <- array(0, c(p, p, df))
    
    centralMatrices[,,1] <- matrices[,,1] - Sigma
    corelMatrics[,,1] <- matrices[,,1] - Sigma
    pelet <- matrices[,,1]
    
    for(i in 2:df){
      centralMatrices[,,i] <- matrices[,,i] - Sigma # epsilon[i]
      corelMatrics[,,i] <- centralMatrices[,,i] # X[i] = epsilon[i]
      
      if(!is.null(AR)){
        arlag <- min(maxAR, i - 1)
        corelMatrics[,,i] <- corelMatrics[,,i] + summatrix(corelMatrics, (i - arlag):(i - 1), rev(AR[1:arlag]))
        # X[i] = epsilon[i] + X[i-1] + ...
      }
      if(!is.null(MA)){
        malag <- min(maxMA, i - 1)
        corelMatrics[,,i] <- corelMatrics[,,i] + summatrix(centralMatrices, (i - malag):(i - 1) , rev(MA[1:malag]))
        # X[i] = epsilon[i] + X[i-1] + ... + epsilon[i-1] + ...
      }
      
      pelet <- pelet + corelMatrics[,,i] + Sigma
      # X[i] = epsilon[i] + X[i-1] + ... + epsilon[i-1] + ... + mu
      # W = sum(X[i])
    }
    
    return(pelet)
  }
  
  return(simplify2array(mclapply(1:n, rawFun, mc.cores = ncores)))
}
