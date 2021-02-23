cov_known_mu <- function(x, y = NULL, mu, mu_x, mu_y, na.rm = FALSE, use = "everything"){

  if(is.vector(x)){
    if(na.rm) x <- x[!is.na(x)]
    if(is.null(y)) {
      y <- x
      if(missing(mu_x)) mu_x <- mu
      mu_y <- mu_x
    }
    out <- mean((x - mu_x) * (y - mu_y))
  }
  else if (is.matrix(x) | is.data.frame(x)){
    stop('var_known_mu supports univariate') # currently univariate
  } else {
    stop('x must be a vector, matrix a data.frame')
  }
  
  return(out)
}

var_known_mu <- function(x, mu = 0, na.rm = FALSE, use = "everything")
  cov_known_mu(x = x, y = NULL, mu_x = mu, na.rm = na.rm, use = use)

sd_known_mu <- function(x, mu, na.rm = FALSE)
  sqrt(var_known_mu(if (is.vector(x) || is.factor(x)) x else as.double(x), mu = mu, na.rm = na.rm))
