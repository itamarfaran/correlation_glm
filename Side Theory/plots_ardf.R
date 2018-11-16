library(tidyverse)

fx <- function(Ti, phi, Limit = TRUE) Ti - (1 + phi)/(1 - phi) + (!Limit)*2*phi*(1 - phi^Ti)/(Ti*(1 - phi)^2)

data.frame(phi = seq(-1, 1, length.out = 1000)) %>% mutate(DF = fx(1, phi, TRUE)) %>%
  ggplot(aes(x = phi, y = DF)) +  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_line(col = 'blue', size = 1) +
  ylim(-10, 2) + labs(title = "DF as a function of Phi, T = 1", x = "Phi", y = "Effect on DF")

data.frame(phi = seq(0, 1, length.out = 1000)) %>% mutate(Real = fx(100, phi, FALSE), Limit = fx(100, phi, TRUE)) %>%
  gather(key = 'Type', value = 'DF', -phi) %>% ggplot(aes(x = phi, y = DF, col = Type)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + geom_line(size = 1) + ylim(-150, 150)
  labs(title = "DF as a function of Phi, T = 1", x = "Phi", y = "Effect on DF")


f2x <- function(Ti, phi) -2*phi*(1 - phi^Ti)/((1 - phi)^2 * Ti)

data.frame(Ti = seq(1, 250, by = 1)) %>%
  mutate('-0.75' = f2x(Ti, -0.75), '-0.5' = f2x(Ti, -0.5),
         '-0.25' = f2x(Ti, -0.25), '0' = f2x(Ti, 0), '0.25' = f2x(Ti, 0.25), '0.5' = f2x(Ti, 0.5),
         '0.75' = f2x(Ti, 0.75)) %>%
  gather(key = "Phi", value = "Value", -Ti) %>% ggplot(aes(x = Ti, y = Value, col = Phi)) + ylim(-1, 1) +
  geom_line(size = 0.65) + labs(title = "Residual as a function of T and Phi", x = "T", y = "Residual")

f3x <- function(Ti, phi) Ti - (1 + phi)/(1 - phi) + 2*phi*(1 - phi^Ti)/(Ti*(1 - phi)^2)

create_table <- function(maxTi, FUN, name){
  data.frame(Ti = seq(1, maxTi, 1)) %>%
    mutate('0.1' = FUN(Ti, 0.1), '0.5' = FUN(Ti, 0.5), '0.8' = FUN(Ti, 0.8),
           '0.95' = FUN(Ti, 0.95)) %>% gather(key = "Phi", value = "DF", -Ti) %>% mutate(Type = as.character(name))
}

rbind( create_table(1000, fx, 'Limit'), create_table(1000, f3x, 'Real') ) %>% 
  ggplot(aes(x = Ti, y = DF, col = Phi, linetype = Type)) + geom_hline(yintercept = 0) + geom_line() + scale_x_continuous(trans = 'log10')

rbind( create_table(150, fx, 'Limit'), create_table(150, f3x, 'Real') ) %>% 
  ggplot(aes(x = Ti, y = DF, col = Phi, linetype = Type)) + geom_hline(yintercept = 0) + geom_line()

rbind( create_table(50, fx, 'Limit'), create_table(50, f3x, 'Real') ) %>% 
  ggplot(aes(x = Ti, y = DF, col = Phi, linetype = Type)) + geom_hline(yintercept = 0) + geom_line()
