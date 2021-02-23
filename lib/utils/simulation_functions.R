# todo: need methodology for simulating alphas with dim>1
# todo: why not package as well?

build_parameters <- function(p, percent_alpha, range_alpha, dim_alpha = 1, loc_scale = c(0,1), enforce_min_alpha=FALSE, seed = NULL){
  #Build Real Sigma and Theta
  if(!is.null(seed)) set.seed(seed[1])
  tmp_theta <- rnorm(2*p^2, loc_scale[1], loc_scale[2])
  tmp_theta <- matrix(tmp_theta, nrow = 2*p)
  tmp_theta <- t(tmp_theta) %*% tmp_theta
  real_theta <- force_symmetry(cov2cor(tmp_theta))
  
  #Generate Variances and create Variance matrix parameters
  if(!is.null(seed)) set.seed(seed[2])
  sds <- abs(rlogis(p))
  real_sigma <- diag(sds) %*% real_theta %*% diag(sds)
  
  #Build Real Alpha
  sum_alpha <- rep(1, p)
  n_alpha <- floor(percent_alpha * p)
  non_null_alphas <- sample(p, n_alpha)
  
  sum_alpha[non_null_alphas] <- runif(n_alpha, range_alpha[1], range_alpha[2])
  if(enforce_min_alpha) sum_alpha[sample(non_null_alphas, 1)] <- range_alpha[1]
  
  alpha <- matrix(runif(p*dim_alpha), nr = p)
  alpha <- apply(alpha, 1, function(x) x/sum(x))
  
  if(dim_alpha > 1) alpha <- t(alpha)
  alpha <- alpha * (sum_alpha %o% rep(1, dim_alpha))
  
  return(list(corr_mat = real_theta, cov_mat = real_sigma, alpha = alpha))
}


create_samples <- function(n_sim = 1, n_h, n_s, p, Tlength = 100,
                           percent_alpha = 0, range_alpha = 0:1, enforce_min_alpha=FALSE, dim_alpha = 1, loc_scale = c(0,1),
                           linkFun = linkFunctions$multiplicative_identity, real_theta = NULL, real_sick = NULL,
                           ARsick = NULL, ARhealth = NULL, MAsick = NULL, MAhealth = NULL, random_effect = NULL,
                           seed = NULL, ncores = 1){
  
  parameters <- build_parameters(
    p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, dim_alpha = dim_alpha,
    loc_scale = loc_scale, enforce_min_alpha = enforce_min_alpha, seed = seed
    )
  if(is.null(real_theta)) real_theta <- parameters$corr_mat
  if(is.null(real_sick)){
    alpha <- parameters$alpha
    g11 <- linkFun$FUN(t = triangle2vector(real_theta), a = alpha, d = dim_alpha)
  } else {
    alpha <- NULL
    g11 <- real_sick
  }
  
  rawFun <- function(b){
    list(
      healthy = create_correlation_matrices(
        real_corr = real_theta, sample_size = n_h, df = Tlength,
        AR = ARhealth, MA = MAhealth, random_effect = random_effect, ncores = 1
        ),
      sick = create_correlation_matrices(
        real_corr = g11, sample_size = n_s, df = Tlength,
        AR = ARsick, MA = MAsick, random_effect = random_effect, ncores = 1
        )
    )
  }
  
  output <- list(real_theta = real_theta, alpha = alpha, samples = list())
  
  output$samples <- if(n_sim == 1) rawFun(1) else mclapply(1:n_sim, rawFun, mc.cores = ncores)
  
  return(output)
}

