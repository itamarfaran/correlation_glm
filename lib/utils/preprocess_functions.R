calculate_na_mat <- function(array_){
  na_mat <- matrix(0, nr=dim(array_)[1], nc=dim(array_)[2])
  for(i in seq_len(dim(array_)[3])){
    na_mat <- na_mat + is.na(array_[,,i])
  }
  na_mat <- na_mat/i
  return(na_mat)
}


get_percent_na <- function(na_mat){
  return(na_mat[1,])
}


which_na_greater_threshold <- function(na_mat, threshold){
  which(get_percent_na(na_mat) >= threshold)
}


get_count_na_by_threshold <- function(na_mat){
  possible_thresholds <- sort(unique(get_percent_na(na_mat)))
  percent_over_threshold <- numeric(length(possible_thresholds))
  for(i in seq_along(possible_thresholds))
    percent_over_threshold[i] <- length(which_na_greater_threshold(na_mat, possible_thresholds[i]))
  percent_over_threshold <- percent_over_threshold/length(get_percent_na(na_mat))
  return(cbind(possible_thresholds, percent_over_threshold))
}


filter_columns_by_na_threshold <- function(array_, threshold){
  which_columns_keep <- which(get_percent_na(calculate_na_mat(array_)) <= threshold)
  which_columns_to_drop <- which(get_percent_na(calculate_na_mat(array_)) > threshold)
  print(which_columns_to_drop)  # todo: should be part of output!
  new_array_ <- array_[which_columns_keep, which_columns_keep, ]
  return(new_array_)
}


who_to_drop <- function(array_){
  return(unique(which(is.na(array_), arr.ind = T)[,3]))
}


drop_na_by_percent <- function(array_, threshold=NULL, verbose=FALSE){
  if(is.null(threshold)){
    na_mat <- calculate_na_mat(array_)
    threshold_tbl <- get_count_na_by_threshold(na_mat)
    colnames(threshold_tbl) <- c('threshold', 'percent_columns_ommited')
    
    threshold_vect <- threshold_tbl[,1]
    threshold <- min(threshold_vect[threshold_vect > 0])
    if(verbose) print(data.frame(threshold_tbl))
    message(paste0('taking ', round(threshold, 5), ' as threshold'))
  }
  orig_dim <- dim(array_)
  array_new_ <- filter_columns_by_na_threshold(array_, threshold)
  message(paste0('orig p = ', dim(array_)[1], ', return p = ', dim(array_new_)[1]))
  return(array_new_)
}


prepare_corrmat_data <- function(link, corr_matrix_name, healthy_index_name, sick_index_name, subset=NULL, threshold=NULL){
  real_dta <- R.matlab::readMat(link)
  corr_mats <- real_dta[[corr_matrix_name]]
  
  for(try in 1:10){
    if(class(corr_mats) == 'array') break()
    corr_mats <- simplify2array(corr_mats)
    if(try == 10) stop('corr_mats could not be simplified to array')
  }
  
  data_and_list <- drop_na_by_percent(corr_mats, threshold=threshold)
  reduced_data <- data_and_list
  which_cols_na <- unique(which(is.na(reduced_data), T)[,3])
  
  for(i in seq_len(dim(reduced_data)[3])) reduced_data[,,i] <- force_symmetry(reduced_data[,,i])
  
  healthy_index <- as.vector(real_dta[[healthy_index_name]])
  healthy_index <- healthy_index[!(healthy_index %in% which_cols_na)]
  sick_index <- as.vector(real_dta[[sick_index_name]])
  sick_index <- sick_index[!(sick_index %in% which_cols_na)]
  
  healthy_dta <- reduced_data[,,healthy_index]
  sick_dta <- reduced_data[,,sick_index]
  
  if(!is.null(subset)){
    healthy_dta <- healthy_dta[subset, subset, ]
    sick_dta <- sick_dta[subset, subset, ]
  }
  
  p <- dim(sick_dta)[1]
    
  sample_data <- list(samples = list(healthy = healthy_dta, sick = sick_dta),
                      p = p, which_cols_na = which_cols_na, link = link)
  return(sample_data)
}

