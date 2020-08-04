##### source and load #####
source("main_work/simulations/auxilary_functions.R")
ipak('ggExtra', 'plotly')

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

transpose(multiplicative_load$estimates$steps)$alpha %>% 
  lapply(as.vector) %>% do.call(cbind, .) %>%
  data.table() -> alpha_steps

colnames(alpha_steps) <- as.character(1:ncol(alpha_steps))

alpha_steps[,index := 1:.N]
alpha_steps <- melt(alpha_steps, id.vars = 'index')
alpha_steps[,variable := as.numeric(variable)]
alpha_steps[,index := factor(index)]

samp <- alpha_steps[,sample(index, 20)]
p <- ggplot(alpha_steps[index %in% samp], aes(x = variable, y = value, col = index)) +
  geom_hline(yintercept = 1) + geom_line() + geom_point() +
  theme(legend.position = 'none'); p

ggplotly(p)

transpose(quotent_load$estimates$steps)$value %>% 
  do.call(c, .) %>%
  data.table(value = .) %>%
  (function(x) {x[,index:=1:.N]; x}) %>%
  (function(x) x[index>1]) %>%
  ggplot(aes(x = index, y = value)) + 
  geom_point()
