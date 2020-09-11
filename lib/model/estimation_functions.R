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

corrmat_covariance <- function(matr, nonpositive = c("stop", "force", "ignore"), fisher_z = FALSE, use_cpp = FALSE){
  if(is.vector(matr)) matr <- vector2triangle(matr, diag_value = 1)
  
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
    cppFunction(paste0(scan(
        "lib/model/corcalc_c.cpp",
        what = "character",
        sep = "\n",
        quiet = TRUE),
        collapse = "\n"))

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
  
  if(fisher_z){
    diagonalized <- triangle2vector(matr)
    transformed <- 1/(1 - diagonalized^2)
    gradient <- diag(transformed)
    output <- gradient %*% output %*% gradient
  }
  return(output)
}

corrmat_covariance_from_dt <- function(dt, use_cpp = FALSE, est_n = FALSE, only_diag = TRUE, ncores = 1){
  out <- mclapply(
    seq_len(nrow(dt)),
    function(i) corrmat_covariance(dt[i,], 'ignore', use_cpp = use_cpp),
    mc.cores = ncores
    )
  
  sigma <- calculate_mean_matrix(simplify2array(out))
  
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


inner_optim_loop <- function(
  healthy_dt, sick_dt, alpha0 = NULL, theta0 = NULL, weight_matrix = NULL,
  dim_alpha = 1, linkFun = linkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iter_config = list(), optim_config = list(), early_stop = TRUE,
  verbose = TRUE){
  
  if('reltol' %in% names(iter_config) & 'abstol' %in% names(iter_config))
    stop('can supply only one of reltol or abstol')
  
  model_reg_config <- modifyList(list(lambda = 0, lp = 2), model_reg_config)
  matrix_reg_config <- modifyList(list(method = 'constant', const = 0), matrix_reg_config)
  iter_config <- modifyList(list(max_loop = 50, reltol = 1e-06, min_loop = 3), iter_config)
  optim_config <- modifyList(list(method = "BFGS", reltol = 1e-06, log_optim = FALSE), optim_config)

  p <- 0.5 + sqrt(1 + 8*ncol(sick_dt))/2
  m <- 0.5*p*(p-1)
  
  if(is.null(theta0)) theta0 <- colMeans(healthy_dt)
  if(is.null(alpha0)) alpha0 <- matrix(linkFun$NULL_VAL, nr = p, nc = dim_alpha)

  if(!is.positive.semi.definite(vector2triangle(theta0, diag_value = 1))
    || !is.positive.semi.definite(linkFun$FUN(t = theta0, a = alpha0, d = dim_alpha)))
    warning("Initial parameters dont result with positive-definite matrices")

  dim_alpha <- length(alpha0)/p
  if(dim_alpha %% 1 != 0) stop("alpha0 not multiplicative of p")

  if(is.null(weight_matrix)) {
    weight_matrix <- weight_matrix_reg <- weight_matrix_reg_inv <- diag(m)
  } else {
    weight_matrix_reg <- regularize_matrix(
      weight_matrix,
      method = matrix_reg_config$method,
      const = matrix_reg_config$const)
      weight_matrix_reg_inv <- solve(weight_matrix_reg)
  }

  temp_theta <- theta0
  temp_alpha <- alpha0
  log_optim_out <- list()
  steps <- list()
  convergence <- NA_integer_
  steps[[1]] <- list(
    theta = temp_theta,
    alpha = temp_alpha,
    value = sum_of_squares(
      theta = temp_theta,
      alpha = temp_alpha,
      sick_dt = sick_dt,
      inv_sigma = weight_matrix_reg_inv,
      linkFun = linkFun,
      dim_alpha = dim_alpha,
      reg_lambda = model_reg_config$lambda,
      reg_p = model_reg_config$lp,
    )
  )


  tt <- Sys.time()
  if(verbose) message(paste0("Time of intialization: ", tt, "; Progress: 'Loop, (Time, Convergence, Distance)'"))
  for(i in 2:iter_config$max_loop){

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
      inv_sigma = weight_matrix_reg_inv,
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
    
    temp_alpha <- optim_alpha$par
    steps[[i]] <- list(
      theta = temp_theta,
      alpha = temp_alpha,
      value = optim_alpha$value,
      convergence = optim_alpha$convergence
      )
    convergence <- c(convergence, optim_alpha$convergence)
    log_optim_out[[i]] <- if(optim_config$log_optim) optim_alpha else NA

    # Stopping rule
    if('abstol' %in% names(iter_config)){
      distance <- sqrt(mean((steps[[i]]$alpha - steps[[i-1]]$alpha)^2))
      distance_lower_than_threshold <-  distance < iter_config$abstol
    } else {
      distance <- abs(steps[[i-1]]$value - steps[[i]]$value)
      distance_lower_than_threshold <- distance < (iter_config$reltol * (abs(steps[[i]]$value) + iter_config$reltol))
    }

    if(verbose) cat(paste0(
      i, " (", round(as.double.difftime(Sys.time() - tt, units = "secs")), "s, ",
      steps[[i]]$convergence, ", ", round(distance, 5), "); "
    ))

    condition0 <- FALSE
    if(i > iter_config$min_loop){
      look_back <- iter_config$min_loop - 1
      index <- if(look_back > 0) i - 0:look_back else i
      condition0 <- distance_lower_than_threshold & (sum(convergence[index]) == 0)
    }
    
    # Early Stop
    if(early_stop){
      did_converge <- convergence[length(convergence)] == 0
      did_minimize <- steps[[i]]$value <= steps[[i-1]]$value
      
      if(did_converge & !did_minimize){
        steps[[i]] <- NULL
        i <- i - 1
        
        temp_theta <- steps[[i]]$theta
        temp_alpha <- steps[[i]]$alpha
        convergence[length(convergence)] <- -1
        condition0 <- TRUE
        warning('early stopping used; last iteration didn\'t minimize target')
      }
    }
      
    if(condition0) break()
  }
  if(i == iter_config$max_loop) warning('optimization reached maximum iterations')
  if(verbose){
    tt <- Sys.time() - tt
    units(tt) <- "secs"
    tt <- as.numeric(tt)
    message(paste0("\nTotal time: ", floor(tt/60), " minutes and ", round(tt %% 60, 1), " seconds."))
  }

  output <- list(
    theta = temp_theta,
    alpha = temp_alpha,
    linkFun = linkFun,
    vcov = weight_matrix_reg_inv,
    convergence = convergence,
    steps = steps, log_optim = log_optim_out
    )
  
  return(output)
}


estimate_alpha <- function(
  healthy_dt, sick_dt, dim_alpha = 1,
  linkFun = linkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  raw_start = TRUE, iid_config = list(), cov_config = list(),
  bias_correction = TRUE, early_stop = TRUE, verbose = TRUE){
  
  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }
  if(length(verbose) == 1) verbose <- rep(verbose, 2)
  
  healthy_dt <- convert_corr_array_to_data_matrix_test(healthy_dt)
  sick_dt <- convert_corr_array_to_data_matrix_test(sick_dt)
  
  alpha0 <- theta0 <- NULL
  if(!raw_start){
    iid_model <- inner_optim_loop(
      healthy_dt = healthy_dt, sick_dt = sick_dt,
      alpha0 = NULL, theta0 = NULL,
      weight_matrix = NULL, dim_alpha = dim_alpha,
      linkFun = linkFun,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = iid_config$iter_config, optim_config = iid_config$optim_config,
      early_stop = early_stop, verbose = verbose[1]
    )
    alpha0 <- iid_model$alpha
    theta0 <- iid_model$theta
  }

  weight_matrix <- corrmat_covariance_from_dt(sick_dt)

  cov_model <- inner_optim_loop(
    healthy_dt = healthy_dt, sick_dt = sick_dt,
    alpha0 = alpha0, theta0 = theta0,
    weight_matrix = weight_matrix, linkFun = linkFun,
    model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
    iter_config = cov_config$iter_config, optim_config = cov_config$optim_config,
    early_stop = early_stop, verbose = verbose[2]
  )
  
  if(bias_correction) cov_model$alpha <- cov_model$alpha - median(cov_model$alpha) + linkFun$NULL_VAL
  
  return(cov_model)
}


estimate_alpha_jacknife <- function(
  healthy_dt, sick_dt, dim_alpha = 1, 
  alpha0 = NULL, theta0 = NULL,
  linkFun = linkFunctions$multiplicative_identity,
  model_reg_config = list(), matrix_reg_config = list(),
  iid_config = list(iter_config = list(min_loop = 0)), cov_config = list(),
  return_gee = FALSE, jack_healthy = TRUE, early_stop = TRUE,
  verbose = TRUE, ncores = 1
  ){
  
  mclapply_ <- if(verbose) pbmcapply::pbmclapply else parallel::mclapply
  
  apply_fun <- function(i, boot_dt){
    if(boot_dt == 'sick'){
      sick_dt_ <- sick_dt[-i,]
      healthy_dt_ <- healthy_dt
    } else if(boot_dt == 'healthy'){
      sick_dt_ <- sick_dt
      healthy_dt_ <- healthy_dt[-i,]
    }
    
    weight_matrix <- corrmat_covariance_from_dt(sick_dt_)
    
    cov_model <- inner_optim_loop(
      healthy_dt = healthy_dt_, sick_dt = sick_dt_,
      alpha0 = alpha0, theta0 = theta0,
      weight_matrix = weight_matrix, linkFun = linkFun,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = cov_config$iter_config, optim_config = cov_config$optim_config,
      early_stop = early_stop, verbose = FALSE
    )

    gee_out <- if(return_gee){
      triangle2vector(
        compute_gee_variance(
          cov_obj = cov_model,
          healthy_dt = healthy_dt_,
          sick_dt = sick_dt_,
          est_mu = TRUE
        ),
        diag = TRUE
      )
    } else NA
    
    return(list(
      theta = cov_model$theta,
      alpha = cov_model$alpha,
      convergence = tail(cov_model$convergence, 1), 
      gee_var = gee_out
    ))
  }
  
  for(name in c('iter_config', 'optim_config')){
    if(!name %in% names(iid_config)) iid_config[[name]] <- list()
    if(!name %in% names(cov_config)) cov_config[[name]] <- list()
  }
  
  healthy_dt <- convert_corr_array_to_data_matrix_test(healthy_dt)
  sick_dt <- convert_corr_array_to_data_matrix_test(sick_dt)
  
  if(is.null(alpha0) | is.null(theta0)){
    iid_model <- inner_optim_loop(
      healthy_dt = healthy_dt, sick_dt = sick_dt,
      alpha0 = alpha0, theta0 = theta0,
      weight_matrix = NULL, dim_alpha = dim_alpha,
      linkFun = linkFun,
      model_reg_config = model_reg_config, matrix_reg_config = matrix_reg_config,
      iter_config = iid_config$iter_config, optim_config = iid_config$optim_config,
      early_stop = early_stop, verbose = FALSE
    )
    if(is.null(alpha0)) alpha0 <- iid_model$alpha
    if(is.null(theta0)) theta0 <- iid_model$theta
  }
  
  if(verbose) cat('\njacknifing Sick Observations...\n')
  cov_obj_sick <- mclapply_(seq_len(nrow(sick_dt)), apply_fun, boot_dt = 'sick', mc.cores = ncores)
  cov_obj_sick_t <- purrr::transpose(cov_obj_sick)
  
  theta <- do.call(rbind, cov_obj_sick_t$theta)
  alpha <- do.call(rbind, lapply(cov_obj_sick_t$alpha, as.vector))
  convergence <- do.call(c, cov_obj_sick_t$convergence)
  gee_var <- do.call(rbind, cov_obj_sick_t$gee_var)

  if(jack_healthy){
    if(verbose) cat('\njacknifing Healthy Observations...\n')
    cov_obj_healthy <- mclapply_(seq_len(nrow(healthy_dt)), apply_fun, boot_dt = 'healthy', mc.cores = ncores)
    cov_obj_healthy_t <- purrr::transpose(cov_obj_healthy)
    
    theta_h <- do.call(rbind, cov_obj_healthy_t$theta)
    alpha_h <- do.call(rbind, lapply(cov_obj_healthy_t$alpha, as.vector))
    convergence_h <- do.call(c, cov_obj_healthy_t$convergence)
    gee_var_h <- do.call(rbind, cov_obj_healthy_t$gee_var)
    
    theta <- rbind(theta_h, theta)
    alpha <- rbind(alpha_h, alpha)
    convergence <- c(convergence_h, convergence)
    gee_var <- rbind(gee_var, gee_var_h)
    is_sick <- c(rep(0, nrow(healthy_dt)), rep(1, nrow(sick_dt)))
  }
  
  return(list(
    theta = theta,
    alpha = alpha,
    gee_var = gee_var,
    convergence = convergence,
    is_sick = is_sick,
    linkFun = linkFun
  ))
}
