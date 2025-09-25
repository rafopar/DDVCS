
void DrawAccEffects(){

  TCanvas *c1 = new TCanvas("c1", "", 950, 950);
  
  TFile *file_in = new TFile("AnaGrape_22GeV_18.root", "Read");

  TH2D *h_Q2_xB2 = (TH2D*)file_in->Get("h_Q2_xB2");
  h_Q2_xB2->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
  
  const double Eb = 22.;
  const double Mp = 0.9383;
  const double th_min = 7.5;
  const double P_min = 1.;

  TF1 *f_constTh = new TF1("f_constTh", "2*[0]*[1]*x/( 1 + [0]*x/([1]*(1-cos([2]*TMath::DegToRad() ))) )", 0., 1. );
  f_constTh->SetParameters(Mp, Eb, th_min);
  f_constTh->SetLineWidth(3);

  TF1 *f_constP = new TF1("f_constP", "2*[0]*([1] - [2])*x", 0., 1);
  f_constP->SetLineColor(6);
  f_constP->SetParameters(Mp, Eb, P_min);
  f_constP->SetLineWidth(3);
  
  h_Q2_xB2->Draw();
  f_constTh->Draw("Same");
  f_constP->Draw("Same");

  TLatex *lat1 = new TLatex();
  lat1->SetNDC();
  lat1->SetTextFont(42);
  lat1->SetTextColor(2);
  lat1->DrawLatex(0.45, 0.5, "#theta = 7.5 deg");
  lat1->SetTextColor(6);
  lat1->SetTextAngle(75);
  lat1->DrawLatex(0.2, 0.65, "P = 1 GeV");
  c1->Print("Figs/Q2_xB_Acc_Effects.pdf");
  c1->Print("Figs/Q2_xB_Acc_Effects.png");
  c1->Print("Figs/Q2_xB_Acc_Effects.root");
  
}
