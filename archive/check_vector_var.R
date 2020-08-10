p <- 7

A <- (vector_to_triangle(1:(p*(p-1)/2)) - diag(p))

vector_var_matrix_calc_COV <- function(MATR){
  p <- dim(MATR)[1]
  m <- (p*(p-1)/2)
  
  v1 <- numeric(0)
  v2 <- numeric(0)
  for(i in 1:(p-1)){
    v1 <- c(v1, rep(i, p-i))
    v2 <- c(v2, (i+1):p)
  }
  order_vect <- cbind(v1,v2)
  
  pelet <- matrix(nrow = m, ncol = m)
  for(i in 1:m){
    for(j in i:m){
      indexes <- c(order_vect[i,], order_vect[j,])
      pelet[i,j] <- as.integer(paste0(indexes[1],indexes[2],indexes[3],indexes[4]))
      pelet[j,i] <- pelet[i,j]
    }
  }
  
  return(pelet)
}

B1 <- vector_var_matrix_calc_COV(A)
B1
A %>% triangle_to_vector()

vector_var_matrix_calc_COR <- function(MATR, nonpositive = c("Stop", "Force", "Ignore"),
                                       reg_par = 0){
  if(length(nonpositive) > 1) nonpositive <- nonpositive[1]
  if(!is.positive.definite(MATR)){
    if(nonpositive == "Force") {MATR <- force_positive_definiteness(MATR)$Matrix
    } else if(nonpositive != "Ignore") stop("MATR not positive definite") }
  
  real.cov2 <- function(i, j, k, l) {
    (MATR[i,j]*MATR[k,l]/2) * (MATR[i,k]^2 + MATR[i,l]^2 + MATR[j,k]^2 + MATR[j,l]^2) -
      MATR[i,j]*(MATR[i,k]*MATR[i,l] + MATR[j,k]*MATR[j,l]) -
      MATR[k,l]*(MATR[i,k]*MATR[j,k] + MATR[i,l]*MATR[j,l]) +
      (MATR[i,k]*MATR[j,l] + MATR[i,l]*MATR[j,k])
  }
  
  p <- dim(MATR)[1]
  m <- (p*(p-1)/2)
  
  v1 <- numeric(0)
  v2 <- numeric(0)
  for(i in 1:(p-1)){
    v1 <- c(v1, rep(i, p-i))
    v2 <- c(v2, (i+1):p)
  }
  order_vect <- cbind(v1,v2)
  
  i <- 3
  j <- 7
  
  pelet <- matrix(nrow = m, ncol = m)
  for(i in 1:m){
    for(j in i:m){
      indexes <- c(order_vect[i,], order_vect[j,])
      #pelet[i,j] <- real.cov2(indexes[1], indexes[2], indexes[3], indexes[4])
      pelet[i,j] <- as.integer(paste0(indexes[1],indexes[2],indexes[3],indexes[4]))
      pelet[j,i] <- pelet[i,j]
    }
  }
  
  if((reg_par<0)|(reg_par>1)) warning("Regularization Parameter not between 0,1")
  pelet <- (1 - reg_par)*pelet + reg_par*diag(diag(pelet))
  
  return(pelet)
}

B2 <- vector_var_matrix_calc_COR(A, nonpositive = "Ignore")
B2
A %>% triangle_to_vector()
