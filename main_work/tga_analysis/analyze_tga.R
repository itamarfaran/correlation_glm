##### source and load #####
source("main_work/simulations/auxilary_functions.R")
ipak('ggExtra')

w <- 300
h <- 300

load('main_work/tga_analysis/analysis_tga_multiplicative_identity.RData')
multiplicative_load <- list(
  estimates = results,
  variance = gee_var
)
load('main_work/tga_analysis/analysis_tga_additive_quotent.RData')
quotent_load <- list(
  estimates = results,
  variance = gee_var
)



##### plot explanatory #####
png('main_work/tga_analysis/control_explanatory.png', w, h)
corrplot(
  corr = calculate_mean_matrix(sample_data$samples$healthy), 
  method = 'color', tl.pos = "n")
dev.off()

png('main_work/tga_analysis/control_explanatory_diag.png', 2*w, 2*h)
corrplot(
  corr = calculate_mean_matrix(sample_data$samples$healthy), 
  method = 'color', tl.pos = "n", type = 'upper', diag = F)
dev.off()

png('main_work/tga_analysis/diagnosed_explanatory.png', w, h)
corrplot(
  corr = calculate_mean_matrix(sample_data$samples$sick), 
  method = 'color', tl.pos = "n")
dev.off()

col_pal <- colorRampPalette(c('#E3C000', '#B04500', '#000000', '#003672', '#00DD0A'))

empirical_difference_corrmat <- with(sample_data$samples, calculate_mean_matrix(sick) - calculate_mean_matrix(healthy))
png('main_work/tga_analysis/difference_explanatory.png', w, h)
corrplot(
  corr = empirical_difference_corrmat, 
  method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35))
dev.off()

png('main_work/tga_analysis/difference_explanatory_diag.png', 2*w, 2*h)
corrplot(
  corr = empirical_difference_corrmat, 
  method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35), type = 'upper', diag = F)
dev.off()


##### analyze t-tests #####
fisher_z <- function(r) 0.5*log((1 + r)/(1 - r))
t_test_results <- matrix(0, sample_data$p, sample_data$p)
for(i in 1:(sample_data$p - 1)) for(j in (i + 1):sample_data$p) t_test_results[i,j] <- with(
  sample_data$samples, t.test(fisher_z(healthy[i,j,]), fisher_z(sick[i,j,]))$p.value)
t_test_results[upper.tri(t_test_results)] <- p.adjust(t_test_results[upper.tri(t_test_results)], 'BH')
t_test_results <- t_test_results + t(t_test_results)
diag(t_test_results) <- 1
print(range(t_test_results))


##### organize gee results #####
end_results <- data.table(
  index = 1:sample_data$p,
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


to_write <- end_results_long[,.(index, method, est = round(est, 2), p_adjusted = round(p_adjusted, 3))]
to_write[,p_adjusted := case_when(
  p_adjusted == 0 ~ '(<.001)',
  TRUE ~ paste0('(', p_adjusted, ')')
)]
to_write[,sig_codes := case_when(
  p_adjusted > .1 ~ '   ',
  p_adjusted > .05 ~ ' . ',
  p_adjusted > .01 ~ ' * ',
  p_adjusted > .001 ~ ' **',
  TRUE ~ '***'
)]
to_write[,`:=`(
  Index = index,
  method = case_when(
    method == 'multiplicative' ~ 'Multiplicative',
    method == 'quotent' ~ 'Quotent'
  ),
  Estimate = paste0(est, ' ', sig_codes),
  Pvalue = p_adjusted,
  index = NULL, est = NULL, p_adjusted = NULL, sig_codes = NULL
  )]

to_write <- melt(to_write, id.vars = c('Index', 'method'), value.name = 'Estimate')[order(Index, method)]
to_write <- dcast(to_write, Index + variable  ~ method)
to_write[,`:=`(variable = NULL, Index = as.character(Index), temp = 1:.N)]
to_write[temp %% 2 == 0, Index := '']
to_write[,temp:=NULL]
fwrite(to_write, file = 'main_work/tga_analysis/results.csv')



##### plot estimates #####
estimates_plt <- dcast(end_results_long[,.(index, method, est)], index ~ method) %>% 
  ggplot(aes(x = multiplicative, y = quotent)) + 
  geom_vline(xintercept = 1, linetype = 1, size = 1, col = 'darkgrey') + 
  geom_hline(yintercept = 0, linetype = 1, size = 1, col = 'darkgrey') + 
  geom_point(alpha = 0.8, shape = 21, color = 'black', fill = 'grey') + 
  labs(
    title = 'Estimates under Different Link Functions',
    x = TeX('$\\Lambda_{d,ij}=\\Theta_{ij}\\alpha_{i}\\alpha_{j}$'),
    y = TeX('$\\Lambda_{d,ij}=\\Theta_{ij}/\\left(1+\\alpha_{i}+\\alpha_{j}\\right)$')
  ) + 
  theme_user()


zscores_plt <- dcast(end_results_long[,.(index, method, z_value)], index ~ method) %>% 
  ggplot(aes(x = multiplicative, y = quotent)) + 
  geom_vline(xintercept = 0, linetype = 1, size = 1, col = 'darkgrey') + 
  geom_hline(yintercept = 0, linetype = 1, size = 1, col = 'darkgrey') + 
  geom_point(alpha = 0.8, shape = 21, color = 'black', fill = 'grey') + 
  xlim(-5, 5) + ylim(-5, 5) + 
  labs(
    title = 'Z-Scores under Different Link Functions',
    x = TeX('$\\Lambda_{d,ij}=\\Theta_{ij}\\alpha_{i}\\alpha_{j}$'),
    y = TeX('$\\Lambda_{d,ij}=\\Theta_{ij}/\\left(1+\\alpha_{i}+\\alpha_{j}\\right)$')
  ) + 
  theme_user()

zscores_plt <- ggMarginal(zscores_plt, type = "histogram", fill = 'grey', col = 'white')

custom_ggsave('estimates_links.png', estimates_plt)
custom_ggsave('zscores_links.png', zscores_plt)

mod <- quotent_load
estimate_difference_corrmat_quot <- with(mod$estimates, linkFun$FUN(theta, alpha, 1) -
                                           vector2triangle(theta, diag_value = 1))

png('main_work/tga_analysis/difference_model.png', w, h)
corrplot(
  corr = estimate_difference_corrmat_quot, 
  method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.15, .15))
  #cl.lim = c(-.35, .35))
dev.off()


alpha_effect <- with(mod$estimates, linkFun$FUN(rep(1, length(theta)), alpha, 1))
alpha_effect_pmat <- end_results_long[method == 'quotent', p_adjusted] %o% rep(1, sample_data$p)
for(i in 1:sample_data$p) for (j in i:sample_data$p)
  alpha_effect_pmat[i,j] <- alpha_effect_pmat[j,i] <- min(alpha_effect_pmat[i,j], alpha_effect_pmat[j,i])

png('main_work/tga_analysis/alpha_effect_all.png', w, h)
corrplot(corr = alpha_effect - 1, 
         method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35),
         p.mat = alpha_effect_pmat, insig = "n", bg = 'lightgrey')
dev.off()
png('main_work/tga_analysis/alpha_sig.png', w, h)
corrplot(corr = alpha_effect - 1, 
         method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35),
         p.mat = alpha_effect_pmat, insig = "blank", bg = 'lightgrey')
dev.off()



# diagnosed_diffence <- calculate_mean_matrix(sample_data$samples$sick) - 
#   with(mod$estimates, linkFun$FUN(theta, alpha, 1))
# 
# plot(with(mod$estimates, linkFun$FUN(theta, alpha, 1)), diagnosed_diffence)
# 
# corrplot(
#   corr = diagnosed_diffence, 
#   method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100))#, cl.lim = c(-.15, .15))


