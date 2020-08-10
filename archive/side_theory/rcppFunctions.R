library(Rcpp)



pelet1 <- vector_var_matrix_calc_COR(MATR)

cppFunction('NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
            NumericVector out(nrow);
            
            for (int i = 0; i < nrow; i++) {
            double total = 0;
            for (int j = 0; j < ncol; j++) {
            total += x(i, j);
            }
            out[i] = total;
            }
            return out;
            }')

cppFunction('
bool allC(NumericVector x){
  int n = x.size();
  for(int i = 0; i < n; ++i){
    if(!x[i]) return false;
  }
  return true;
}
')

x <- as.logical(sample(0:1, 10, T))
allC(x)
all(x)
