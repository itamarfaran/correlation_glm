source("main_work/Code/01_generalFunctions.R")
source("main_work/Code/02_simulationFunctions.R")
source("main_work/Code/03_estimationFunctions.R")
source("main_work/Code/04_inferenceFunctions.R")

link <- "main_work/Data/NMDA_all_data_AAL90.mat"
linkFun <- linkFunctions$multiplicative_identity
sub_samp_size <- 10:11
boot_size <- 3

Real.dta <- readMat(link)

comparison_function <- function(sub_samp_size, linkFun, seed = NULL){
  which.drop <- which(is.na(Real.dta$group.all[1,,1]), arr.ind = TRUE)
  p <- dim(Real.dta$group.all)[1] - length(which.drop)
  All.data <- array(dim = c(p, p, dim(Real.dta$group.all)[3]))
  
  if(length(which.drop)){
    for(i in 1:dim(All.data)[3])
      All.data[,,i] <- force_symmetry(Real.dta$group.all[-which.drop,-which.drop,i])
  } else {
    for(i in 1:dim(All.data)[3])
      All.data[,,i] <- force_symmetry(Real.dta$group.all[,,i])
  }
  if(!is.null(seed)) set.seed(seed); sub_samp <- sort(sample(p, sub_samp_size))
  
  healthy.Real_t <- All.data[sub_samp, sub_samp, Real.dta$CONTROLS]
  healthy.Real <- if(sum(is.na(healthy.Real_t))) {
    healthy.Real_t[,,-unique(which(is.na(healthy.Real_t), arr.ind = TRUE)[,3])]
  } else healthy.Real_t
  sick.Real_t <- All.data[sub_samp, sub_samp, Real.dta$NMDA]
  sick.Real <- if(sum(is.na(sick.Real_t))) {
    sick.Real_t[,,-unique(which(is.na(sick.Real_t), arr.ind = TRUE)[,3])]
  } else sick.Real_t
  sampleData <- list(healthy = healthy.Real, sick = sick.Real, p = p, which.drop = which.drop, link = link)
  
  all_positive_definite <- 
    abind(sampleData$healthy, sampleData$sick) %>% apply(3, is.positive.definite) %>% all()
  no_na <- sum(is.na(sampleData$healthy)) + sum(is.na(sampleData$sick)) == 0
  
  if(!all_positive_definite | !no_na) stop("stop")
  
  output_ss <- Estimate.Loop(iniAlpha = 0.8, healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                             linkFun = linkFun)
  
  output_ml_full <- 
    estimateAlpha(healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                  T_thresh = 320, linkFun = linkFun, updateU = 1,
                  var_weights = c(1, 0, 0),
                  progress = FALSE)
  
  output_ml_cheat <- 
    estimateAlpha(healthy.data = sampleData$healthy, sick.data = sampleData$sick,
                  T_thresh = 320, linkFun = linkFun, updateU = 1,
                  var_weights = c(0, 1, 0),
                  progress = FALSE)
  
  results <- list(output_ss, output_ml_cheat, output_ml_full)
  
  likelihoods <- 
    sapply(results,
           function(out) with(out,
                              minusloglik(
                                theta = theta, alpha = alpha,
                                sick.data = cor.matrix_to_norm.matrix(sampleData$sick))
           )
    )
  
  alphas <- sapply(results, function(l) l$alpha)
  
  names(likelihoods) <- c("sum_squares", "lik_partial", "lik_full") 
  colnames(alphas) <- c("sum_squares", "lik_partial", "lik_full") 
  
  return(list(likelihoods = likelihoods,
              alphas = alphas))
}
res <- list()
alpha_list <- list()
alpha_dist_list <- list()
likelihood_list <- list()

for(i in 1:length(sub_samp_size)){
  sss <- sub_samp_size[i]
  current_models <- mclapply(1:boot_size,
                               function(b) comparison_function(sss, linkFun),
                               mc.cores = ncores)
  alphas <- do.call(rbind, transpose(current_models)$alphas)
  alphas <- cbind(rep(sss, sss*boot_size), alphas)
  names(alphas)[1] <- "dim"
  
  alphas_distance <- 
    t(sapply(transpose(current_models)$alphas, function(mat) c(
      vnorm(mat[,"sum_squares"] - mat[,"lik_partial"], sqrt = TRUE),
      vnorm(mat[,"sum_squares"] - mat[,"lik_full"], sqrt = TRUE)
    )
    ))
  alphas_distance <- cbind(rep(sss, boot_size), alphas_distance)
  colnames(alphas_distance) <- c("dim", "lik_partial", "lik_full")
  
  likelihoods <- do.call(rbind, transpose(current_models)$likelihoods)
  likelihoods <- cbind(rep(sss, boot_size), likelihoods)
  colnames(likelihoods)[1] <- "dim"
  
  res[[i]] <- current_models
  alpha_list[[i]] <- alphas
  alpha_dist_list[[i]] <- alphas_distance
  likelihood_list[[i]] <- likelihoods
}

alpha_list <- as.data.table(do.call(rbind, alpha_list))
alpha_dist_list <- as.data.table(do.call(rbind, alpha_dist_list))
likelihood_list <- as.data.table(do.call(rbind, likelihood_list))

alphas_plot <- 
  alpha_list %>%
    gather(key = model, value = alphas, -dim) %>%
    ggplot(aes(x = factor(dim), y = alphas, fill = model)) +
    geom_boxplot(alpha = 0.1) + 
    geom_dotplot(binaxis = "y", stackdir = "center",
                 # binpositions="all", stackgroups = TRUE,
                 binwidth = 0.01, alpha = 0.7)

alpha_list %>%
  gather(key = model, value = alphas, -sum_squares, -dim) %>% 
  ggplot(aes(x = sum_squares, y = alphas, col = model, shape = factor(dim))) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + 
  geom_abline(slope = 1, intercept = 0)

likelihood_list[,chisq_val := pchisq(2*(lik_partial - lik_full), 1, lower.tail = F)]

likelihood_list %>%
  ggplot(aes(x = chisq_val, fill = factor(dim))) + 
  geom_histogram(binwidth = 1)

save.image("main_work/data/enviroments/comparison_likelihoods.RData")