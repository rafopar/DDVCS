

double calc_Xsec(double charge, double NDet){

  const double rho = 0.07085; // g/cm3
  const double l = 5.; //cm
  const double N_A = 6.02214076e23;
  const double q_e = 1.6e-19; // Electron's charge
  const double pb = 1e-36;
  charge = charge*1e-3; // to conver to miliColumbs
  
  double Lumi = rho*l*N_A*charge/q_e;
  
  double x_sec = (1./pb)*NDet/Lumi; // X-sec in [pb]
  return x_sec; //pb
}
