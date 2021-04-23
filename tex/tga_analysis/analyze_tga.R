##### source and load #####
source("tex/simulations/aux_.R")
ipak('ggExtra')
sig_level = .05

w <- 300
h <- 300

SD_CONST = 1.1

load('tex/tga_analysis/analysis_tga_multiplicative_identity.RData')
multiplicative_load <- list(
  estimates = results,
  variance = gee_var
)
load('tex/tga_analysis/analysis_tga_additive_quotent.RData')
quotent_load <- list(
  estimates = results,
  variance = gee_var
)


##### plot explanatory #####
png('tex/tga_analysis/control_explanatory.png', w, h)
corrplot(
  corr = calculate_mean_matrix(sample_data$samples$healthy), 
  method = 'color', tl.pos = "n")
dev.off()

png('tex/tga_analysis/control_explanatory_diag.png', 2*w, 2*h)
corrplot(
  corr = calculate_mean_matrix(sample_data$samples$healthy), 
  method = 'color', tl.pos = "n", type = 'upper', diag = F)
dev.off()

png('tex/tga_analysis/diagnosed_explanatory.png', w, h)
corrplot(
  corr = calculate_mean_matrix(sample_data$samples$sick), 
  method = 'color', tl.pos = "n")
dev.off()

col_pal <- colorRampPalette(c('#E3C000', '#B04500', '#000000', '#003672', '#00DD0A'))

empirical_difference_corrmat <- with(sample_data$samples, calculate_mean_matrix(sick) - calculate_mean_matrix(healthy))
png('tex/tga_analysis/difference_explanatory.png', w, h)
corrplot(
  corr = empirical_difference_corrmat, 
  method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35))
dev.off()

png('tex/tga_analysis/difference_explanatory_diag.png', 2*w, 2*h)
corrplot(
  corr = empirical_difference_corrmat, 
  method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35), type = 'upper', diag = F)
dev.off()


##### analyze t-tests #####
fisher_z <- function(r) 0.5*log((1 + r)/(1 - r))

t_test_results <- t_test_results_corrected <- matrix(0, sample_data$p, sample_data$p)

for(i in 1:(sample_data$p - 1)) for(j in (i + 1):sample_data$p) t_test_results[i,j] <- with(
  sample_data$samples, t.test(fisher_z(healthy[i,j,]), fisher_z(sick[i,j,]))$p.value)

t_test_results_corrected[upper.tri(t_test_results_corrected)] <- p.adjust(t_test_results[upper.tri(t_test_results)], 'BH')

t_test_results <- t_test_results + t(t_test_results)
t_test_results_corrected <- t_test_results_corrected + t(t_test_results_corrected)

diag(t_test_results_corrected) <- diag(t_test_results) <- 1
print(range(t_test_results))

create_mannahtan_plot <- function(mat, upper_lim = 0){
  mat_dt <- data.table(mat)
  colnames(mat_dt) <- as.character(1:ncol(mat_dt))
  mat_dt[,j := 1:.N]
  mat_dt_long <- melt(mat_dt, id.vars = 'j', variable.name = 'i', value.name = 'p_value')
  mat_dt_long[,i:=as.integer(i)]
  
  out <- ggplot(mat_dt_long, aes(x = i, y = p_value)) + 
    geom_point(alpha = .1)
  if(upper_lim > 0){
    out <- out + scale_y_continuous(trans = reverselog_trans(), limits = c(1, upper_lim)) 
  } else {
    out <- out + scale_y_continuous(trans = reverselog_trans())
  }
    
  return(out)
}
mannahtan_plot <- create_mannahtan_plot(t_test_results_corrected, upper_lim = .1) +
  labs(title = 'P-values From Mass Univariate T-Tests',
       x = '', y = 'P-values (-log10 Scale)') +
  theme_user() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    )
custom_ggsave('mannahtan_plot.png', mannahtan_plot, width = 2)


##### organize gee results #####
add_na_columns <- function(df, indices){
  for(i in indices) df[index >= i, index := index + 1]
  df_out <- rbind(df, data.table(index = indices), fill = TRUE)
  setorder(df_out, index)
  return(df_out)
}

end_results <- add_na_columns(
  data.table(
    index = 1:sample_data$p,
    est_multiplicative = as.vector(multiplicative_load$estimates$alpha),
    sd_multiplicative = sqrt_diag(multiplicative_load$variance)*SD_CONST,
    est_quotent = as.vector(quotent_load$estimates$alpha),
    sd_quotent = sqrt_diag(quotent_load$variance)*SD_CONST
  ),
  c(21, 28, 75, 76)
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
end_results_long[!is.na(est),p_adjusted := p.adjust(p_value, 'BH'), by = method]


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
fwrite(to_write, file = 'tex/tga_analysis/results.csv')



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

mod_name <- c('multiplicative', 'quotent')[1]
mod <- switch(mod_name, 'multiplicative' = multiplicative_load, 'quotent' = quotent_load, NA)

estimate_difference_corrmat_quot <- with(mod$estimates, linkFun$FUN(theta, alpha, 1) -
                                           vector2triangle(theta, diag_value = 1))

png('tex/tga_analysis/difference_model.png', w, h)
corrplot(
  corr = estimate_difference_corrmat_quot, 
  method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35))
  #cl.lim = c(-.35, .35))
dev.off()


alpha_effect <- with(mod$estimates, linkFun$FUN(rep(1, length(theta)), alpha, 1))
alpha_effect_pmat <- end_results_long[method == 'quotent', p_adjusted] %o% rep(1, sample_data$p)
for(i in 1:sample_data$p) for (j in i:sample_data$p)
  alpha_effect_pmat[i,j] <- alpha_effect_pmat[j,i] <- min(alpha_effect_pmat[i,j], alpha_effect_pmat[j,i])

png('tex/tga_analysis/alpha_effect_all.png', w, h)
corrplot(corr = alpha_effect - 1, 
         method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35),
         p.mat = alpha_effect_pmat, insig = "n", bg = 'lightgrey', sig.level = sig_level)
dev.off()

png('tex/tga_analysis/alpha_sig.png', w, h)
corrplot(corr = alpha_effect - 1, 
         method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35),
         p.mat = alpha_effect_pmat, insig = "blank", bg = 'white', sig.level = sig_level)
dev.off()

alpha_effect_pmat_toplot <- alpha_effect_pmat
temp_j <- 1
for(j in seq_len(nrow(alpha_effect_pmat_toplot)))
  if(any(alpha_effect_pmat_toplot[temp_j,] > .05))
    break()

is_sig <- alpha_effect_pmat_toplot[temp_j,] < .05
alpha_effect_pmat_toplot[,!is_sig] <- 1
png('tex/tga_analysis/alpha_sig_var.png', w, h)
corrplot(corr = alpha_effect - 1, 
         method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100), cl.lim = c(-.35, .35),
         p.mat = alpha_effect_pmat_toplot, insig = "blank", bg = 'white', sig.level = sig_level)
dev.off()

find_bh_threshold <- function(pvalues, sig_level = .05){
  pvalues <- sort(pvalues)
  tresholds <- sig_level*seq_along(pvalues)/length(pvalues)
  threshold <- tresholds[max(which(pvalues <= tresholds))]
  return(threshold)
}

find_holmes_threshold <- function(pvalues, sig_level = .05){
  pvalues <- sort(pvalues)
  tresholds <- sig_level/(length(pvalues) + 1 - seq_along(pvalues))
  threshold <- tresholds[min(which(pvalues > tresholds)) - 1]
  return(threshold)
}

find_bonferroni_threshold <- function(pvalues, sig_level = .05){
  return(sig_level/length(pvalues))
}

toplot_manhattan <- end_results_long[method == mod_name]

sum_rejected <- toplot_manhattan[,sum(p_adjusted <= sig_level)]
toplot_manhattan[,`:=`(
  ci_low = est - qnorm(1 - sig_level/sum_rejected/2)*sd*(p_adjusted <= sig_level),
  ci_upp = est + qnorm(1 - sig_level/sum_rejected/2)*sd*(p_adjusted <= sig_level),
  is_significant = ifelse(p_adjusted <= sig_level, 'Significant', 'Insignificant'),
  difference = sign(z_value*(p_adjusted <= sig_level))
)]
toplot_manhattan[,difference := case_when(
  difference == 1 ~ 'Increase',
  difference == -1 ~ 'Decay',
  difference == 0 ~ 'Insignificant',
  TRUE ~ NA_character_
)]

alpha_estimate_plot <- ggplot(toplot_manhattan, aes(x = index, y = est, ymin = ci_low, ymax = ci_upp, shape = difference, fill = difference, size = difference)) + 
  geom_hline(yintercept = mod$estimates$linkFun$null_value, linetype = 2, color = 'darkgrey') + 
  geom_point() + 
  # geom_text(aes(label = index)) + 
  labs(x = '', y = 'Estimate', color = '') +
  theme_user() + 
  theme(legend.position = 'none') +
  scale_shape_manual(values = c('Increase' = 24, 'Decay' = 25, 'Insignificant' = 21)) + 
  scale_fill_manual(values = c('Increase' = 'darkgrey', 'Decay' = 'darkgrey', 'Insignificant' = 'white')) + 
  scale_size_manual(values = c('Increase' = 2.5, 'Decay' = 2.5, 'Insignificant' = 2))


labels_dt <- data.table(
  x = max(toplot_manhattan$index) - 20,
  y = c(
    find_bh_threshold(toplot_manhattan$p_value)*1.6,
    find_holmes_threshold(toplot_manhattan$p_value)*0.6
    )
  )
labels_dt[,label := paste0(c('BH', 'Bonferroni'))]

alpha_manhattan_plot <- ggplot(toplot_manhattan, aes(x = index, y = p_value)) + 
  geom_hline(yintercept = find_bh_threshold(toplot_manhattan$p_value)) + 
  geom_hline(yintercept = find_bonferroni_threshold(toplot_manhattan$p_value)) + 
  geom_point(aes(shape = difference, fill = difference, size = difference)) + 
  scale_y_continuous(
    trans = reverselog_trans(),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  geom_label(aes(x = x, y = y, label = label), data = labels_dt, color = 'black') +
  labs(x = 'Index', y = 'P-value', color = '') +
  theme_user() + 
  theme(legend.position = 'none') + 
  scale_shape_manual(values = c('Increase' = 24, 'Decay' = 25, 'Insignificant' = 21)) + 
  scale_fill_manual(values = c('Increase' = 'darkgrey', 'Decay' = 'darkgrey', 'Insignificant' = 'white')) + 
  scale_size_manual(values = c('Increase' = 2.5, 'Decay' = 2.5, 'Insignificant' = 2))

alpha_manhattan_plot

out <- arrangeGrob(alpha_estimate_plot, alpha_manhattan_plot, nrow = 2)
custom_ggsave('alpha_manhattan.png', out, width = 1.5)

# diagnosed_diffence <- calculate_mean_matrix(sample_data$samples$sick) - 
#   with(mod$estimates, linkFun$FUN(theta, alpha, 1))
# 
# plot(with(mod$estimates, linkFun$FUN(theta, alpha, 1)), diagnosed_diffence)
# 
# corrplot(
#   corr = diagnosed_diffence, 
#   method = 'color', is.corr = F, tl.pos = "n", col = col_pal(100))#, cl.lim = c(-.15, .15))


