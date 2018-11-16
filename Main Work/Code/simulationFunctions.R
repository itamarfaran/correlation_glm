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

create_correlation_matrices <- function(real_corr, sample_size, df = 0, var_scale = NULL,
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
  
  if(length(seed.control)==1) set.seed(seed.control)
  pelet_matrices <- rWishart(sample_size, df, real_var)
  for(b in 1:sample_size) pelet_matrices[,,b] <- force_symmetry(cov2cor(pelet_matrices[,,b]))
  
  return(pelet_matrices)
}
