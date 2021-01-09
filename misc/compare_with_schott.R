H_ <- function(i, j, m, n){
  out <- matrix(0, m, n)
  out[i, j] <- 1
  return(out)
}

K_ <- function(m, n){
  out <- matrix(0, m*n, m*n)
  for(i in 1:m) for(j in 1:n) out <- out + H_(i, j, m, n) %x% t(H_(i, j, m, n))
  return(out)
}

e_ <- function(i, m){
  out <- numeric(m)
  out[i] <- 1
  return(out)
}

E_ <- function(i, j, m){
  out <- e_(i, m) %o% e_(j, m)
  return(out)
}

N_ <- function(m){
  out <- .5 * (diag(m^2) + K_(m, m))
  return(out)
}

LAMBDA_ <- function(m){
  out <- matrix(0, m^2, m^2)
  for(i in 1:m) out <- out + E_(i, i, m) %x% E_(i, i, m)
  return(out)
}

THETA_ <- function(P){
  m <- ncol(P)
  I_m <- diag(m)
  LAMBDA_m <- LAMBDA_(m)
  out <- P %x% P - (I_m %x% P) %*% LAMBDA_m %*% (P %x% P) - 
    (P %x% P) %*% LAMBDA_m %*% (I_m %x% P) + 
    (I_m %x% P) %*% LAMBDA_m %*% (P %x% P) %*% LAMBDA_m %*% (I_m %x% P)
  return(out)
}

var_vec_R <- function(R, n){
  m <- ncol(R)
  out <- N_(m) %*% THETA_(R) %*% N_(m)
  # out <- out * 2 / (n - 1)
  out <- out * 2
  return(out)
}

corrmat_covariance <- function(matr){
  p <- nrow(matr)
  m <- p*(p-1)/2
  order_vecti <- unlist(lapply(1:(p - 1), function(i) rep(i, p - i)))
  order_vectj <- unlist(lapply(1:(p - 1), function(i) (i + 1):p))
  
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
  
  output <- output + t(output) - diag(diag(output))
    return(output)
}


### Testing ###

p <- 90
mat <- matrix(0, p^2, 3)

k <- 1
for(i in 1:p){
  for(j in 1:p){
    mat[k,1] <- k
    mat[k,2] <- i
    mat[k,3] <- j
    k <- k + 1
  }
}
index <- mat[mat[,3] > mat[,2],1]


tt <- matrix(2*runif(p^2) - 1, ncol = p) 
M <- cov2cor(t(tt) %*% tt)
M <- (M + t(M)) / 2

now <- Sys.time()
faran <- corrmat_covariance(M)
faran_runtime <- Sys.time() - now

# now <- Sys.time()
# schott <- var_vec_R(M, 2)
# schott_runtime <- Sys.time() - now
# 
# faran
# schott[index, index]
# 
# max(abs(schott[index, index] - faran))

faran_runtime
# schott_runtime
