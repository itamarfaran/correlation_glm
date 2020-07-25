source("main_work/code/01_general_functions.R")
source("main_work/code/02_simulation_functions.R")
source("main_work/code/03_estimation_functions.R")
source("main_work/code/04_inference_functions.R")
ipak('ggcorrplot', 'GGally')

load('main_work/simulations/analysis_tga_multiplactive.RData')
multiplicative_load <- list(
  estimates = results,
  variance = gee_var
)
load('main_work/simulations/analysis_tga_quotent.RData')
quotent_load <- list(
  estimates = results,
  variance = gee_var
)
raw_data <- sample_data

ggcorr(NULL, cor_matrix = calculate_mean_matrix(raw_data$samples$healthy))

corrplot(
  corr = calculate_mean_matrix(raw_data$samples$healthy), 
  method = 'color', title = 'Control Subjects',
  tl.pos = "n")

corrplot(
  corr = calculate_mean_matrix(raw_data$samples$sick), 
  method = 'color', title = 'Diagnosed Subjects',
  tl.pos = "n")

# todo: different color pallet here
corrplot(
  corr = calculate_mean_matrix(raw_data$samples$sick) - calculate_mean_matrix(raw_data$samples$healthy), 
  method = 'color', title = 'Difference Between Groups', is.corr = F, tl.pos = "n")
  
fisher_z <- function(r) 0.5*log((1 + r)/(1 - r))
t_test_results <- matrix(0, raw_data$p, raw_data$p)
for(i in 1:(raw_data$p - 1)) for(j in (i + 1):raw_data$p) t_test_results[i,j] <- t.test(
    fisher_z(raw_data$samples$healthy[i,j,]),
    fisher_z(raw_data$samples$sick[i,j,])
    )$p.value
t_test_results[upper.tri(t_test_results)] <- p.adjust(t_test_results[upper.tri(t_test_results)], 'BH')
t_test_results <- t_test_results + t(t_test_results)
diag(t_test_results) <- 1

end_results <- data.table(
  index = 1:raw_data$p,
  est_multiplicative = as.vector(multiplicative_load$estimates$alpha),
  sd_multiplicative = sqrt_diag(multiplicative_load$variance),
  est_quotent = as.vector(quotent_load$estimates$alpha),
  sd_quotent = sqrt_diag(quotent_load$variance)
)

end_results_long <- melt(end_results, id.vars = 'index')
cols <- c('type', 'method')
end_results_long[,(cols) := asplit(do.call(rbind, str_split(end_results_long[,variable], '_')), 2)]
end_results_long[,variable := NULL]
end_results_long <- dcast(end_results_long, index + method ~ type)[order(method, index)]
end_results_long[,null_value := ifelse(
  method == 'quotent',
  quotent_load$estimates$linkFun$NULL_VAL,
  multiplicative_load$estimates$linkFun$NULL_VAL
  )]
end_results_long[,z_value := (est - null_value)/sd]
end_results_long[,p_value := 2*pnorm(abs(z_value), lower.tail = F)]
end_results_long[,p_adjusted := p.adjust(p_value, 'BH'), by = method]
