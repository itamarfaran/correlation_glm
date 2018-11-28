build_parameters <- function(p, percent_alpha, range_alpha, loc_scale = c(0,1), seed = NULL){
  #Build Real Sigma and Theta
  if(length(seed) > 0) set.seed(seed[1])
  temp <- matrix(rnorm(2*p^2, loc_scale[1], loc_scale[2]), nrow = 2*p)
  temp <- t(temp)%*%temp
  real.theta <- force_symmetry(cov2cor(temp))
  
  #Generate Variances and create Variance matrix parameters
  if(length(seed) > 0) set.seed(seed[2])
  varss <- rlogis(p)^2
  real.sigma <- sqrt(diag(varss)) %*% real.theta %*% sqrt(diag(varss))
  
  #Build Real Alpha
  alpha <- rep(1,p)
  if(length(seed) > 0) set.seed(seed[3])
  alpha[sample(1:p, floor(percent_alpha * p))] <- runif(floor(percent_alpha * p), range_alpha[1], range_alpha[2])
  
  return(list(Corr.mat = real.theta, Cov.mat = real.sigma, Alpha = alpha))
}

create_correlation_matrices <- function(real_corr, sample_size, df = 0, AR = NULL, MA = NULL, var_scale = NULL,
                                        seed.control = NULL){
  if(!is.positive.definite(real_corr)) stop("real_corr not positive definite")
  p <- nrow(real_corr)
  if((df > 0) & (df <= p + 1)) warning("df <= p + 1, p + 1 used in wishart simulations")
  df <- max(df, p + 1)
  if(length(var_scale)==0){
    var_scale <- runif(p,10,100)
  } else if(length(var_scale)!=p){
    stop("var_scale not in dimension")
  }
  
  Dhalf <- sqrt(diag(var_scale))
  real_var <- Dhalf%*%real_corr%*%Dhalf
  
  if(is.null(AR) & is.null(MA)){
    if(length(seed.control) == 1) set.seed(seed.control)
    pelet_matrices <- rWishart(sample_size, df, real_var)
  } else {
    if(length(seed.control) == 1) stop("seed.control not available for rWishart_ARMA")
    pelet_matrices <- rWishart_ARMA(sample_size, df, real_var, AR = AR, MA = MA)
  }
  for(b in 1:sample_size) pelet_matrices[,,b] <- force_symmetry(cov2cor(pelet_matrices[,,b]))
  
  return(pelet_matrices)
}

rWishart_ARMA <- function(n = 1, df, Sigma, AR = NULL, MA = NULL){
  
  if(!is.null(AR)){
    if(!checkInv(AR)) stop("AR process not stationary.")
    maxAR <- length(AR)
    
  } 
  if(!is.null(MA)){
    if(!checkInv(MA)) stop("MA process not invertable")
    maxMA <- length(MA)
  }
  
  p <- ncol(Sigma)
  blocksDF <- floor(df/p)
  if(blocksDF*p != df) warning(paste0("Using ", blocksDF*p, " Degrees of Freedom instead of ", df))
  
  rawFun <- function(k, blocksDF, p, Sigma, AR = NULL, MA = NULL){
  
    matrices <- rWishart(blocksDF, p, Sigma)
    meanMat <- calculate_mean_matrix(matrices)
    
    centralMatrices <- array(0, c(p, p, blocksDF))
    corelMatrics <- array(0, c(p, p, blocksDF))
    pelet <- matrix(0, nrow = p, ncol = p)
  
    corelMatrics[,,1] <- matrices[,,1] - meanMat
    for(i in 2:blocksDF){
      centralMatrices[,,i] <- matrices[,,i] - meanMat
      corelMatrics[,,i] <- meanMat
      
      if(!is.null(AR)){
        arlag <- min(maxAR, i - 1)
        corelMatrics[,,i] <- corelMatrics[,,i] + summatrix(corelMatrics, (i - arlag):(i - 1), rev(AR[1:arlag]))
      }
      if(!is.null(MA)){
        malag <- min(maxMA, i - 1)
        corelMatrics[,,i] <- corelMatrics[,,i] + summatrix(centralMatrices, (i - malag):(i - 1) , rev(MA[1:malag]))
      }
      
      pelet <- pelet + corelMatrics[,,i]
    }
    
    return(pelet)
  }
  
  sapply(1:n, rawFun, blocksDF = blocksDF, p = p, Sigma = Sigma, AR = AR, MA = MA, simplify = "array")
}

createSamples <- function(B = 1, nH, nS, p, Tlength, percent_alpha, range_alpha, loc_scale = c(0,1), 
                          ARsick = NULL, ARhealth = NULL, MAsick = NULL, MAhealth = NULL, seed = NULL){
  
  parameters <- build_parameters(p, percent_alpha, range_alpha, loc_scale, seed)
  real.theta <- parameters$Corr.mat
  real.sigma <- parameters$Cov.mat
  alpha <- parameters$Alpha
  alpha.mat <- create_alpha_mat(alpha)
  
  if(B == 1) return(list(healthy = create_correlation_matrices(real.theta, nH, Tlength,
                                                               AR = ARhealth, MA = MAhealth),
                           sick = create_correlation_matrices(real.theta*alpha.mat, nS, Tlength,
                                                              AR = ARsick, MA = MAsick),
                         real.theta = real.theta, real.sigma = real.sigma, alpha = alpha))
  
  
  pelet <- list(real.theta = real.theta, real.sigma = real.sigma, alpha = alpha, samples = list())
  
  for(b in 1:B){
    pelet$samples[[b]] <-
      list(healthy = create_correlation_matrices(real.theta, nH, Tlength,
                                                 AR = ARsick, MA = MAsick),
           sick = create_correlation_matrices(real.theta*alpha.mat, nS, Tlength,
                                              AR = ARsick, MA = MAsick))
  }
  return(pelet)
}