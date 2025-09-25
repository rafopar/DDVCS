
double Func_xi(double *x, double *par);

void xi_dependence(){

  const double nueMax = 9.;
  TF1 *f_xi = new TF1("f_xi", Func_xi, 0., 9., 2);

  f_xi->SetParameters(6, 9);
  f_xi->SetNpx(4500);
  f_xi->Draw();

  int n_nue = 10;

  for( int i = 0; i < n_nue; i++ ){
    double nue = nueMax*(1. - double(i)/double(n_nue));
    f_xi->SetParameters(nue, 6);
    f_xi->DrawCopy("Same");
  }
  
}

double Func_xi(double *xx, double *par){

  double Q2 = xx[0];
  double nue = par[0];
  double Qp2 = par[1];

  const double Mp = 0.9383;
  
  double xB = Q2/(2*Mp*nue);

  double xi_prime = xB/(2 - xB);

  double xi = xi_prime*(Q2 + Qp2)/Q2;
  
  return xi;
}
