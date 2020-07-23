### todo: need methodology for simulating alphas

linkFunctions <- list(
  "multiplicative_identity" = list(
    FUN = function(t, a, d) {
      a <- matrix(a, nc = d)
      a_mat <- a %*% t(a)
      diag(a_mat) <- 1
      vector2triangle(t, diag_value = 1)*a_mat
      },
    INV = function(a) a,
    CLEAN = function(dt, a, d){
      a <- matrix(a, nc = d)
      a_mat <- a %*% t(a)
      diag(a_mat) <- 1
      a_vect <- triangle2vector(a_mat)
      return(dt / (rep(1, nrow(dt)) %*% t(a_vect)))
    },
    NULL_VAL = 1
  ),
  "additive_identity" = list(
    FUN = function(t, a, d) {
      a <- matrix(a, nc = d)
      a_mat <- a %*% t(a)
      diag(a_mat) <- 0
      vector2triangle(t, diag_value = 1)*(1 + a_mat)
    },
    INV = function(a) a,
    CLEAN = function(dt, a, d){
      a <- matrix(a, nc = d)
      a_mat <- a %*% t(a)
      diag(a_mat) <- 0
      a_vect <- triangle2vector(a_mat)
      return(dt / (1 + rep(1, nrow(dt)) %*% t(a_vect)))
    },
    NULL_VAL = 0
  ),
  "additive_quotent" = list(
    FUN = function(t, a, d) {
      a <- as.vector(a)
      a_mat_t <- replicate(length(a), a)
      a_mat <- a_mat_t + t(a_mat_t)
      diag(a_mat) <- 0
      vector2triangle(t, diag_value = 1)/(1 + a_mat)
      },
    INV = function(a) a,
    CLEAN = function(dt, a, d) {
      a <- as.vector(a)
      a_mat_t <- replicate(length(a), a)
      a_mat <- a_mat_t + t(a_mat_t)
      diag(a_mat) <- 0
      a_vect <- triangle2vector(a_mat)
      return(dt * (1 + rep(1, nrow(dt)) %*% t(a_vect)))
    },
    NULL_VAL = 0
  )
)


build_parameters <- function(p, percent_alpha, range_alpha, dim_alpha = 1, loc_scale = c(0,1), seed = NULL){
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
  sum_alpha[sample(p, n_alpha)] <- runif(n_alpha, range_alpha[1], range_alpha[2])
  
  alpha <- matrix(runif(p*dim_alpha), nr = p)
  alpha <- apply(alpha, 1, function(x) x/sum(x))
  
  if(dim_alpha > 1) alpha <- t(alpha)
  alpha <- alpha * (sum_alpha %o% rep(1, dim_alpha))
  
  return(list(corr_mat = real_theta, cov_mat = real_sigma, alpha = alpha))
}


create_correlation_matrices <- function(real_corr, real_var, sample_size, df = 0, AR = NULL, MA = NULL,
                                        random_effect = NULL, seed_control, ncores = 1){
  if(missing(real_var)){
    if(!is.positive.definite(real_corr)) stop("real_corr not positive definite")
    p <- nrow(real_corr)
    var_scale <- runif(p, 10, 100)
    dhalf <- sqrt(diag(var_scale))
    real_var <- force_symmetry(dhalf %*% real_corr %*% dhalf)
  } else {
    if(!missing(real_corr)) warning("Both real_corr and real_var provided; Using real_var")
  }
  if(!is.positive.definite(real_var)) stop("real_var not positive definite")
  if(is.null(AR) & is.null(MA)){
    if(!missing(seed_control)){
      if(df >= p) {
        set.seed(seed_control)
      } else {
        stop("seed_control not available for rWishart2")
      }
    }
    out <- rWishart2(sample_size, df, real_var, random_effect = random_effect)
  } else {
    if(!missing(seed_control)) stop("seed_control not available for rWishart_ARMA")
    out <- rWishart_ARMA(sample_size, df, real_var, AR = AR, MA = MA, random_effect = random_effect, ncores = ncores)
  }

  return(
    simplify2array(
      lapply(
        1:sample_size,
        function(b)
          force_symmetry(cov2cor(out[,,b]))
        )
      )
    )
}


gen_random_effect_sigma(Sigma, random_effect = NULL){
  if(is.null(random_effect)) return(Sigma)
  if(0 <= random_effect & random_effect < Inf) p <- p*(1 + 1/random_effect)
  return(rWishart(1, p, Sigma)[,,1]/p)
}


rWishart2 <- function(n = 1, df, Sigma, random_effect = NULL){
  p <- ncol(Sigma)
  
  if(df >= p & is.null(random_effect)) return(rWishart(n, df, Sigma))
  if(df < p) warning("Wishart degrees of freedom lower than matrix dimension.")
  
  rawFun <- function(k){
    Sigma <- gen_random_effect_sigma(Sigma, random_effect)
    matrices <- rmvnorm(n = df, sigma = Sigma)
    return(t(matrices) %*% matrices)
  }
  simplify2array(lapply(1:n, rawFun))
}


rWishart_ARMA <- function(n = 1, df, Sigma, AR = NULL, MA = NULL, random_effect = NULL, ncores = 1){
  p <- ncol(Sigma)
  
  if(is.null(MA) & is.null(AR)) return(rWishart2(n = n, df = df, Sigma = Sigma, random_effect = random_effect))
  if(df < p) warning("Wishart degrees of freedom lower than matrix dimension.")
  
  max_ar <- max_ma <- NULL
  
  if(!is.null(AR)){
    if(!check_invertability_arma(AR)) stop("AR process not stationary.")
    max_ar <- length(AR)
  }
  if(!is.null(MA)){
    if(!check_invertability_arma(MA)) stop("MA process not invertable.")
    max_ma <- length(MA)
  }

  rawFun <- function(k){
    Sigma <- gen_random_effect_sigma(Sigma, random_effect)
    normal_matrix <- rmvnorm(n = df, sigma = Sigma)
    normal_matrix_arma <- normal_matrix
    
    for(i in 2:df){
      if(!is.null(AR)){
        arlag <- min(max_ar, i - 1)
        normal_matrix_arma[i,] <- normal_matrix_arma[i,] + vector_sum(normal_matrix_arma, (i - arlag):(i - 1), rev(AR[1:arlag]))
        # X[i] = epsilon[i] + X[i-1] + ...
      }
      if(!is.null(MA)){
        malag <- min(max_ma, i - 1)
        normal_matrix_arma[i,] <- normal_matrix_arma[i,] + vector_sum(normal_matrix, (i - malag):(i - 1) , rev(MA[1:malag]))
        # X[i] = epsilon[i] + X[i-1] + ... + epsilon[i-1] + ...
      }
    }
    return(t(normal_matrix_arma) %*% normal_matrix_arma)
  }
  
  return(simplify2array(mclapply(1:n, rawFun, mc.cores = ncores)))
}


create_samples <- function(n_sim = 1, n_h, n_s, p, Tlength,
                           percent_alpha, range_alpha, dim_alpha = 1, loc_scale = c(0,1),
                           linkFun = linkFunctions$multiplicative_identity,
                           ARsick = NULL, ARhealth = NULL, MAsick = NULL, MAhealth = NULL, random_effect = NULL,
                           seed = NULL, ncores = 1){
  
  parameters <- build_parameters(
    p = p, percent_alpha = percent_alpha, range_alpha = range_alpha, dim_alpha = dim_alpha,
    loc_scale = loc_scale, seed = seed
    )
  
  real_theta <- parameters$corr_mat
  real_sigma <- parameters$cov_mat
  alpha <- parameters$alpha
  g11 <- linkFun$FUN(t = triangle2vector(real_theta), a = alpha, d = dim_alpha)
  
  rawFun <- function(b){
    list(
      healthy = create_correlation_matrices(
        real_corr = real_theta, sample_size = n_h, df = Tlength,
        AR = ARsick, MA = MAsick, random_effect = random_effect, ncores = 1
        ),
      sick = create_correlation_matrices(
        real_corr = g11, sample_size = n_s, df = Tlength,
        AR = ARsick, MA = MAsick, random_effect = random_effect, ncores = 1
        )
    )
  }
  
  output <- list(real_theta = real_theta, real_sigma = real_sigma, alpha = alpha, samples = list())
  
  output$samples <- if(n_sim == 1) list(rawFun(1)) else mclapply(1:n_sim, rawFun, mc.cores = ncores)
  
  return(output)
}

