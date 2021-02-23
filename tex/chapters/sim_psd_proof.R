#### funcs ####

ipak <- function(..., only_install = FALSE){
  pkg <- unlist(list(...))
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  if(!only_install) sapply(pkg, require, character.only = TRUE)
}


force_symmetry <- function(matr)
  return((matr + t(matr))/2)


triangle2vector <- function(matr, diag = FALSE){
  if(nrow(matr) != ncol(matr)) stop("Matrix not p x p")
  return(as.vector(matr[lower.tri(matr, diag = diag)]))
}


vector2triangle <- function(vect, diag = FALSE, diag_value = NA){
  m <- length(vect)
  
  one <- ifelse(diag, -1, 1)
  p <- 0.5*c(one + sqrt(1 + 8*m), one - sqrt(1 + 8*m))
  p <- p[which( (p == round(p)) & p == abs(p) )]
  if(length(p) == 0) stop("Vect length does not fit size of triangular matrix")
  
  out <- matrix(0, ncol = p, nrow = p)
  out[lower.tri(out, diag = diag)] <- vect
  
  out <- out + t(out) - diag(diag(out))
  if(!diag) diag(out) <- diag_value
  return(out)
}


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


#### config ####
ipak('data.table', 'pbmcapply', 'RSpectra', 'matrixcalc')

ncores <- ifelse(tolower(.Platform$OS.type) == "windows", 1, max(detectCores() - 2, 1))

p <- 10
percent_alpha <- 1
range_alpha <- c(-2,2)
n_sim <- 10000

linkfunc <- function(t, a, d) {
  a <- matrix(a, nc = d)
  a_mat <- a %*% t(a)
  diag(a_mat) <- 1
  vector2triangle(t, diag_value = 1)*a_mat
}


#### lapply func ####

simulate <- function(i=NULL, p, percent_alpha, range_alpha, linkfunc){
  params <- build_parameters(p, percent_alpha, range_alpha)
  theta <- force_symmetry(params$corr_mat)
  alpha <- as.vector(params$alpha)
  
  min_eig_val <- eigs_sym(theta, 1, "SM")$values
  alpha_range <- (max(alpha^2) - 1)/min(alpha^2)
  
  g <- linkfunc(triangle2vector(theta), alpha, 1)
  
  is_psd <- is.positive.semi.definite(g)
  is_condition <- min_eig_val >= alpha_range
  
  output <- data.table(sim_id=i,
                      max_alpha=max(alpha),
                      min_alpha=min(alpha),
                      alpha_range=alpha_range,
                      min_eig_val=min_eig_val,
                      is_psd=is_psd,
                      is_condition=is_condition)
  return(output)
}

sim_res <- do.call(rbind, pbmclapply(seq_len(n_sim),
                                     simulate,
                                     p=p,
                                     percent_alpha=percent_alpha,
                                     range_alpha=range_alpha,
                                     linkfunc=linkfunc,
                                     mc.cores=ncores))


#### results ####

print(sim_res)
sim_res[is_condition == TRUE, mean(is_psd)]
sim_res[is_psd == FALSE, mean(is_condition)]

sim_res[is_condition == FALSE, mean(is_psd)]
sim_res[is_psd == TRUE, mean(is_condition)]
