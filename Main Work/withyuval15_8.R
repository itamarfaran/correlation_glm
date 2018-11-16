
Pelet_Cov$Steps[[9]]$alpha

compute_estimated_N_2(, ,
                      , 400)

x <- cov(cor.matrix_to_norm.matrix(sick)) %>% triangle_to_vector(diag = TRUE)
y <- (vector_var_matrix_calc_COR(vector_to_triangle(Pelet_Cov$Steps[[9]]$theta)
                           *create_alpha_mat(Pelet_Cov$Steps[[9]]$alpha))) %>% triangle_to_vector(diag = TRUE)

x <- cov(cor.matrix_to_norm.matrix(sick)) %>% diag()
y <- (vector_var_matrix_calc_COR(vector_to_triangle(Pelet_Cov$Steps[[9]]$theta)
                                 *create_alpha_mat(Pelet_Cov$Steps[[9]]$alpha))) %>% diag()


plot(x,y, pch = 16, col = rgb(0,0.5,1,0.2))

model1 <- lm(y ~ 0 + x)
abline(h = 0)
abline(v = 0)
abline(model1, col = "red")
abline(a = 0, b = 90, col = "blue")

summary(model1)
lm(x~y)
