NumericMatrix corcalc_c(NumericMatrix MATR,
                        int p, int m,
                        NumericVector order_vecti,
                        NumericVector order_vectj)
{
  NumericMatrix pelet(m, m); 

  for (int i1 = 0; i1 < m; i1++) {
    for (int j1 = i1; j1 < m; j1++) {
      int i = order_vecti[i1];
      int j = order_vectj[i1];
      int k = order_vecti[j1];
      int l = order_vectj[j1];
      
      double MATRij = MATR(i,j);
      double MATRkl = MATR(k,l);
      double MATRik = MATR(i,k);
      double MATRil = MATR(i,l);
      double MATRjk = MATR(j,k);
      double MATRjl = MATR(j,l);
      
      pelet(i1,j1) =
        (MATRij*MATRkl/2) * (pow(MATRik, 2) + pow(MATRil, 2) + pow(MATRjk, 2) + pow(MATRjl, 2)) -
        MATRij*(MATRik*MATRil + MATRjk*MATRjl) -
        MATRkl*(MATRik*MATRjk + MATRil*MATRjl) +
        (MATRik*MATRjl + MATRil*MATRjk);
    }
  }
  return(pelet);
}