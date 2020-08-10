---
title: "Univariate Delta Method"
author: "Itamar Faran"
date: "6 April 2018"
output: pdf_document
---






```r
#Builds Sigma
Bui_Sigma <- function(Sx, Sy, rho)  matrix(c(Sx^2, rho*Sx*Sy, rho*Sx*Sy, Sy^2), ncol = 2)

#Creates B samples of size n with covariances
create_samples <- function(B, n, Sx, Sy, rho){
  BIG <- mvrnorm(n = B*n, mu = c(0,0), Sigma = Bui_Sigma(0.5,2,0.5))
  
  pelet <- array(dim = c(n, 2, B))
  for(b in 1:B) pelet[,,b] <- BIG[( (b-1)*n + 1 ):( b*n ), ]
  
  return(pelet)
}

#Raw delta method, eg. compute matrices by computer
delta_meth <- function(n, Sx, Sy, rho){
  Var_mat <- matrix(0, ncol = 3, nrow = 3)
  Var_mat[1,2] <- rho * Sx^3 * Sy
  Var_mat[1,3] <- rho * Sx * Sy^3
  Var_mat[2,3] <- (rho * Sx * Sy)^2
  
  Var_mat <- Var_mat + t(Var_mat)
  Var_mat[1,1] <- 0.5*(1+rho^2) * Sx^2 * Sy^2
  Var_mat[2,2] <- Sx^4
  Var_mat[3,3] <- Sy^4
  Var_mat <- 2*Var_mat/n
  
  grad <- matrix(0, nrow = 3)
  grad[1,1] <- 1/(Sx*Sy)
  grad[2,1] <- -rho/(2*Sx^2)
  grad[3,1] <- -rho/(2*Sy^2)
  
  return(as.vector(t(grad)%*%Var_mat%*%grad))
}

#Delta method we found (final answer only)
delta_meth_simple <- function(n, rho) (1-rho^2)^2/n

print_percent <- function(index, Total, percent) if(((100/percent)*(index/Total))%%1 == 0) cat(100*(index/Total), "%, ")
```

# Simulation

```r
RUN <- TRUE
if(RUN){
  n <- 100
  B <- 5000
  Sx <- sqrt(14)
  Sy <- sqrt(3541)
  rho <- 0.5
  
  B2 <- 1000
  corrs_vars <- numeric(B2)
  for(j in 1:B2){
    Samples <- create_samples(B, n, Sx, Sy, rho)
    corrs_vars[j] <- var(sapply(1:B, function(b) cor(Samples[,1,b], Samples[,2,b])))
  
    print_percent(j, B2, 5)
  }
  
  last_data <- list(n = n, B = B, B2 = B2, Sx = Sx, Sy = Sy, rho = rho, corrs_vars = corrs_vars)
  #save(last_data, file = "last_data")
  rm(last_data)
} else {
  load(file = "last_data")
  
  n <- last_data$n
  B <- last_data$B
  Sx <- last_data$Sx
  Sy <- last_data$Sy
  rho <- last_data$rho
  B2 <- last_data$B2
  corrs_vars <- last_data$corrs_vars

  rm(last_data)
}
```

```
## 5 %, 10 %, 15 %, 20 %, 25 %, 30 %, 35 %, 40 %, 45 %, 50 %, 55 %, 60 %, 65 %, 70 %, 75 %, 80 %, 85 %, 90 %, 95 %, 100 %,
```

# Results and Plots

```r
mean(corrs_vars)
```

```
## [1] 0.005761594
```

```r
delta_meth(n, Sx, Sy, rho)
```

```
## [1] 0.005625
```

```r
delta_meth_simple(n, rho)
```

```
## [1] 0.005625
```

Bias/Real = $\frac{\widehat{Var}\left(\hat{\rho}_{n}\right)-Var\left(\hat{\rho}_{n}\right)}{Var\left(\hat{\rho}_{n}\right)}$

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -0.04246  0.01115  0.02418  0.02428  0.03791  0.10460
```

![](delta_method1_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> ![](delta_method1_files/figure-latex/unnamed-chunk-6-2.pdf)<!-- --> 
