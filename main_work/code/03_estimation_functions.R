compute_estimated_n_raw <- function(est, theo, only_diag = FALSE){
  if(only_diag){
    x <- diag(theo)
    y <- diag(est)
  } else {
    x <- triangle2vector(theo, diag = TRUE)
    y <- triangle2vector(est, diag = TRUE)
  }
  
  return(lm(x ~ 0 + y)$coef)
}

compute_estimated_n <- function(dt, only_diag = TRUE){
  dt <- convert_corr_array_to_data_matrix_test(dt)
  est <- t(dt) %*% dt / nrow(dt)
  theo <- corrmat_covariance_from_dt(dt)
  return(compute_estimated_n_raw(est = est, theo = theo, only_diag = only_diag))
}

corrmat_covariance <- function(matr, nonpositive = c("stop", "force", "ignore"), use_cpp = FALSE){
  
  nonpositive <- match.arg(nonpositive, c("stop", "force", "ignore"))
  if(!is.positive.definite(matr)){
    if(nonpositive == "force") {
      matr <- as.matrix(Matrix::nearPD(matr, corr = TRUE, doSym = TRUE)$mat)
    }
    else if(nonpositive != "ignore") {
      stop("matr not positive definite")
    }
  }
  
  p <- nrow(matr)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i))) - use_cpp
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p)) - use_cpp
  
  if(use_cpp){
    cppFunction(
      paste0(scan(
        "main_work/Code/corcalc_c.cpp",
        what = "character",
        sep = "\n",
        quiet = TRUE),
        collapse = "\n")
    )
    output <- corcalc_c(matr, p, m, order_vecti, order_vectj)
  } else {
    output <- matrix(0, nrow = m, ncol = m)
    for(i1 in 1:m){
      for(j1 in i1:m){
        i <- order_vecti[i1]
        j <- order_vectj[i1]
        k <- order_vecti[j1]
        l <- order_vectj[j1]
        
        matr_ij <- matr[i,j]
        matr_kl <- matr[k,l]
        matr_ik <- matr[i,k]
        matr_il <- matr[i,l]
        matr_jk <- matr[j,k]
        matr_jl <- matr[j,l]
        
        output[i1,j1] <-
          (matr_ij*matr_kl/2) * (matr_ik^2 + matr_il^2 + matr_jk^2 + matr_jl^2) -
          matr_ij*(matr_ik*matr_il + matr_jk*matr_jl) -
          matr_kl*(matr_ik*matr_jk + matr_il*matr_jl) +
          (matr_ik*matr_jl + matr_il*matr_jk)
      }
    }
  }
  
  output <- output + t(output) - diag(diag(output))
  return(output)
}

corrmat_covariance_from_dt <- function(dt, est_n = FALSE, only_diag = TRUE){
  sigma <- corrmat_covariance(vector2triangle(colMeans(dt), diag_value = 1))
  if(est_n){
    est <- t(dt) %*% dt / nrow(dt)
    estimated_n <- compute_estimated_n_raw(est = est, theo = sigma, only_diag = only_diag)
    sigma <- sigma/estimated_n
  }
  return(sigma)
}

theta_of_alpha <- function(alpha, healthy_dt, sick_dt, linkFun, d = 1){
  colMeans(rbind(
    linkFun$CLEAN(dt = sick_dt, a = alpha, d = d),
    healthy_dt
  ))
}


sum_of_squares <- function(
  alpha, theta, sick_dt, inv_sigma,
  linkFun, sigma, dim_alpha = 1,
  reg_lambda = 0, reg_p = 2){
  
  if(missing(inv_sigma)) inv_sigma <- solve(sigma)
  
  g11 <- as.matrix(triangle2vector(linkFun$FUN(t = theta, a = alpha, d = dim_alpha)))
  sse <- nrow(sick_dt) * t(g11) %*% inv_sigma %*% ( 0.5 * g11 - colMeans(sick_dt) )
  
  if(reg_lambda > 0) sse <- sse + reg_lambda*sum((alpha - linkFun$NULL_VAL)^reg_p)
  return(sse)
}


estimate_loop <- function(
  healthy_dt, sick_dt, alpha0 = NULL, theta0 = NULL, dim_alpha = 1,
  linkFun = linkFunctions$multiplicative_identity,
  cov_method = c('identity', 'corrmat'),
  model_reg_config = list(), matrix_reg_config = list(),
  iter_config = list(), optim_config = list(),
  verbose = TRUE){
  
  if('reltol' %in% names(iter_config) & 'abstol' %in% names(iter_config))
    stop('can supply only one of reltol or abstol')
  
  model_reg_config <- modifyList(list(lambda = 0, lp = 2), model_reg_config)
  matrix_reg_config <- modifyList(list(do_reg = FALSE, method = 'constant', const = 1), matrix_reg_config)
  iter_config <- modifyList(list(max_loop = 50, reltol = 1e-06, min_reps = 3), iter_config)
  optim_config <- modifyList(list(method = "BFGS", reltol = 1e-06, log_optim = FALSE), optim_config)
  cov_method <- match.arg(cov_method, c('identity', 'corrmat'))
  
  healthy_dt <- convert_corr_array_to_data_matrix_test(healthy_dt)
  sick_dt <- convert_corr_array_to_data_matrix_test(sick_dt)
  
  healthy_n <- nrow(healthy_dt)
  sick_n <- nrow(sick_dt)
  
  p <- 0.5 + sqrt(1 + 8*ncol(sick_dt))/2
  m <- 0.5*p*(p-1)
  
  if(is.null(theta0)) theta0 <- colMeans(healthy_dt)
  if(is.null(alpha0)) alpha0 <- matrix(1, nr = p, nc = dim_alpha)
  dim_alpha <- length(alpha0)/p
  if(dim_alpha %% 1 != 0) stop("alpha0 not multiplicative of p")
  
  if(!(
    is.positive.definite(vector2triangle(theta0, diag_value = 1)) &
    is.positive.definite(linkFun$FUN(t = theta0, a = alpha0, d = dim_alpha))
  )) warning("Initial parameters dont result with positive-definite matrices")
    
  g12 <- switch(
    cov_method,
    'identity' = diag(m),
    'corrmat' = corrmat_covariance_from_dt(sick_dt),
    NA
    )
  
  g12_reg <- if(matrix_reg_config$do_reg) {
    if(cov_method == 'identity') stop('cannot set do_reg=TRUE with cov_method=\'identity\'')
    regularize_matrix(
      g12,
      method = matrix_reg_config$method,
      const = matrix_reg_config$const
      )
    } else g12
  
  solve_g12_reg <- solve(g12_reg)
  
  temp_theta <- theta0
  temp_alpha <- alpha0
  steps <- list()
  steps[[1]] <- list(theta = temp_theta, alpha = temp_alpha, value = NA)
  log_optim_out <- list()

  #Convergence is a matrix wich tells us if the convergence in each iteration is completed
  convergence <- rep(-1, iter_config$max_loop)
  convergence[1] <- 0
  condition0 <- FALSE
  
  tt <- Sys.time()
  if(verbose) message(paste0("Time of intialization: ", tt, "; Progress: 'Loop, (Time, Convergence, Distance)'"))
  for(i in 2:iter_config$max_loop){
    g11 <- linkFun$FUN(t = temp_theta, a = temp_alpha, d = dim_alpha)
    
    temp_theta <- theta_of_alpha(
      alpha = temp_alpha,
      healthy_dt = healthy_dt,
      sick_dt = sick_dt,
      linkFun = linkFun,
      d = dim_alpha) 
    
    optim_alpha <- optim(
      par = temp_alpha,
      fn = sum_of_squares,
      theta = temp_theta,
      sick_dt = sick_dt,
      inv_sigma = solve_g12_reg,
      linkFun = linkFun,
      dim_alpha = dim_alpha,
      reg_lambda = model_reg_config$lambda,
      reg_p = model_reg_config$lp,
      method = optim_config$method,
      control = list(
        maxit = min(max(500, i*100), 2000),
        reltol = optim_config$reltol
        )
      )
    
    convergence[i] <- optim_alpha$convergence
    temp_alpha <- optim_alpha$par
    steps[[i]] <- list(
      theta = temp_theta,
      alpha = temp_alpha,
      value = optim_alpha$value,
      convergence = convergence[i]
      )
    log_optim_out[[i]] <- if(optim_config$log_optim) optim_alpha else NA
    
    # Stopping rule
    if('abstol' %in% names(iter_config)){
      distance <- sqrt(mean((steps[[i]]$alpha - steps[[i-1]]$alpha)^2))
      distance_lower_than_threshold <-
        distance < iter_config$abstol
    } else {
      distance <- abs(steps[[i-1]]$value - steps[[i]]$value)
      distance_lower_than_threshold <-
        distance < (iter_config$reltol * (abs(steps[[i]]$value) + iter_config$reltol))
    }
    
    if(verbose) cat(paste0(
      i, " (", round(as.double.difftime(Sys.time() - tt, units = "secs")), "s, ",
      convergence[i], ", ", round(distance, 5), "); "
      ))
    
    condition0 <- FALSE
    if(i > iter_config$min_reps) condition0 <-
      distance_lower_than_threshold & (sum(convergence[i - 0:(iter_config$min_reps - 1)]) == 0)
    if(condition0) break()
  }
  if(i == iter_config$max_loop) warning('optimization reached maximum iterations')
  if(verbose){
    tt <- Sys.time() - tt
    units(tt) <- "secs"
    tt <- as.numeric(tt)
    message(paste0("\nTotal time: ", floor(tt/60), " minutes and ", round(tt %% 60, 1), " seconds."))
  }
  
  suppressWarnings({
    max_convergence <- min(
      (min(which(convergence == -1)) - 1),
      iter_config$max_loop
    )
  })
  
  return( list(
    theta = temp_theta,
    alpha = temp_alpha,
    linkFun = linkFun,
    convergence = convergence[1:max_convergence],
    steps = steps, log_optim = log_optim_out) )
}


estimate_alpha <- function(
  healthy_dt, sick_dt, dim_alpha = 1,
  linkFun = linkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iid_config = list(), cov_config = list(),
  verbose = TRUE){
  
  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }
  
  iid_model <- estimate_loop(
    healthy_dt = healthy_dt, sick_dt = sick_dt, dim_alpha = dim_alpha,
    linkFun = linkFunctions$multiplicative_identity,
    cov_method = 'identity',
    model_reg_config = model_reg_config,
    matrix_reg_config = matrix_reg_config,
    iter_config = iid_config$iter_config,
    optim_config = iid_config$optim_config,
    verbose = verbose
    )
  
  cov_model <- estimate_loop(
    healthy_dt = healthy_dt, sick_dt = sick_dt,
    alpha0 = iid_model$alpha, theta0 = iid_model$theta,
    linkFun = linkFunctions$multiplicative_identity,
    cov_method = 'corrmat',
    model_reg_config = model_reg_config,
    matrix_reg_config = matrix_reg_config,
    iter_config = cov_config$iter_config,
    optim_config = cov_config$optim_config,
    verbose = verbose
  )
  
  return(cov_model)
}

# todo: refactor
estimateAlpha_jacknife <- function(
  healthy.data, sick.data, T_thresh = Inf,
  linkFun = linkFunctions$multiplicative_identity,
  dim_alpha = 1, var_weights = c(1, 0, 0), sigma = 1, reg_lambda = 0, reg_p = 2,
  updateU = 1, jack_healthy = TRUE, progress = TRUE, ncores = 1,
  INIconfig = list(iniAlpha = 0.8, MaxLoop = 500, Persic = 0.001, method = "BFGS"),
  FULconfig = list(
    max.loop = 50, epsIter = 2*10^(-3),
    min_reps = 3, nonpositive = "Ignore",
    method = "BFGS", epsOptim = 10^(-5)
    )
  ){
  
  INIconfig <- modifyList(list(iniAlpha = 0.8, MaxLoop = 500, Persic = 0.001, method = "BFGS"), INIconfig)
  FULconfig <- modifyList(list(max.loop = 50, epsIter = 2*10^(-3), min_reps = 3, nonpositive = "Ignore",
                               method = "Nelder-Mead", epsOptim = 10^(-5)), FULconfig)
  
  if(class(healthy.data) == "array") healthy.data <- cor.matrix_to_norm.matrix(healthy.data)
  if(class(sick.data) == "array") sick.data <- cor.matrix_to_norm.matrix(sick.data)
  
  IID <- Estimate.Loop(
    iniAlpha = linkFun$INV(INIconfig$iniAlpha),
    healthy.data = healthy.data, sick.data = sick.data,
    linkFun = linkFun, reg_lambda = reg_lambda, reg_p = reg_p, dim_alpha = dim_alpha,
    MaxLoop = INIconfig$MaxLoop, Persic = INIconfig$Persic, method = INIconfig$method
    )
  cat('\nJack Knifing Sick Observations...\n')
  COV_sick <- pbmclapply(
    1:nrow(sick.data),
    function(i){
      Estimate.Loop2(
        theta0 = IID$theta, alpha0 = IID$alpha,
        healthy.data = healthy.data, sick.data = sick.data[-i,], T_thresh = T_thresh,
        linkFun = linkFun, var_weights = var_weights,
        sigma = sigma, reg_lambda = reg_lambda, reg_p = reg_p,
        nonpositive = FULconfig$nonpositive, max.loop = FULconfig$max.loop,
        epsIter = FULconfig$epsIter, min_reps = FULconfig$min_reps,
        method = FULconfig$method, epsOptim = FULconfig$epsOptim,
        updateU = updateU, progress = FALSE)
    },
    mc.cores = ncores
  )
  
  COV_sick_transposed <- purrr::transpose(COV_sick)
  
  theta <- do.call(rbind, COV_sick_transposed$theta)
  alpha <- do.call(rbind, COV_sick_transposed$alpha)
  est_n <- do.call(rbind, COV_sick_transposed$Est_N)
  
  if(jack_healthy){
    cat('\nJack Knifing Healthy Observations...\n')
    COV_healthy <- pbmclapply(
      1:nrow(healthy.data),
      function(i){
        Estimate.Loop2(
          theta0 = IID$theta, alpha0 = IID$alpha,
          healthy.data = healthy.data[-i,], sick.data = sick.data, T_thresh = T_thresh,
          linkFun = linkFun, var_weights = var_weights,
          sigma = sigma, reg_lambda = reg_lambda, reg_p = reg_p,
          nonpositive = FULconfig$nonpositive, max.loop = FULconfig$max.loop,
          epsIter = FULconfig$epsIter, min_reps = FULconfig$min_reps,
          method = FULconfig$method, epsOptim = FULconfig$epsOptim,
          updateU = updateU, progress = FALSE)
      },
      mc.cores = ncores
    )
    
    COV_healthy_transposed <- purrr::transpose(COV_healthy)
    
    theta <- rbind(
      theta,
      do.call(rbind, COV_healthy_transposed$theta)
    )
    alpha <- rbind(
      alpha,
      do.call(rbind, COV_healthy_transposed$alpha)
    )
    est_n <- rbind(
      est_n,
      do.call(rbind, COV_healthy_transposed$Est_N)
    )
  }
  
  return(list(
    theta = theta,
    alpha = alpha,
    est_n = est_n,
    linkFun = COV_sick[[1]]$linkFun
  ))
}
