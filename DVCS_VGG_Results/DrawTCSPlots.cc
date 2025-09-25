
void DrawTCSPlots(){

  TCanvas *c1 = new TCanvas("c1", "", 950, 950);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.02);

  TGraph *gr_8GeV = new TGraph("TCS_Eg8_Qp2_1p5_TCSPlusBH_xsec.dat", "%lg %lg");
  gr_8GeV->SetLineColor(2);
  gr_8GeV->SetLineStyle(2);

  TGraph *gr_18GeV = new TGraph("TCS_Eg18_Qp2_1p5_TCSPlusBH_xsec.dat", "%lg %lg");
  gr_18GeV->SetLineColor(4);
  gr_18GeV->SetLineStyle(1);

  TGraph *gr_8GeV_TCSOnly = new TGraph("TCS_Eg8_Qp2_1p5_TCSOnly_xsec.dat", "%lg %lg");
  gr_8GeV_TCSOnly->SetLineColor(2);
  gr_8GeV_TCSOnly->SetLineStyle(2);
  gr_8GeV_TCSOnly->SetLineWidth(3);

  TGraph *gr_18GeV_TCSOnly = new TGraph("TCS_Eg18_Qp2_1p5_TCSOnly_xsec.dat", "%lg %lg");
  gr_18GeV_TCSOnly->SetLineColor(4);
  gr_18GeV_TCSOnly->SetLineStyle(1);
  gr_18GeV_TCSOnly->SetLineWidth(3);

  
  TLegend *leg1 = new TLegend(0.5, 0.65, 0.8, 0.89);
  leg1->SetBorderSize(0);
  leg1->AddEntry(gr_8GeV, "8 GeV: TCS + BH");
  leg1->AddEntry(gr_18GeV, "18 GeV: TCS + BH");
  leg1->AddEntry(gr_18GeV_TCSOnly, "18 GeV: TCS");
  leg1->AddEntry(gr_8GeV_TCSOnly, "8 GeV: TCS");
  
  c1->SetLogy();
  TMultiGraph *mtgr1 = new TMultiGraph();
  mtgr1->Add(gr_18GeV);
  mtgr1->Add(gr_8GeV);
  mtgr1->Add(gr_18GeV_TCSOnly);
  mtgr1->Add(gr_8GeV_TCSOnly);
  mtgr1->Draw("AL");
  leg1->Draw();
  mtgr1->SetTitle("TCS x-sec: Q^{'2} = 1.5 GeV^{2}; -t [GeV^{2}]; d#sigma/dQ^{'2}dt pb/[GeV^{4}]");
  c1->Print("Figs/TCS_Qp2_1p5.pdf");
  c1->Print("Figs/TCS_Qp2_1p5.png");
  c1->Print("Figs/TCS_Qp2_1p5.root");
  
}
