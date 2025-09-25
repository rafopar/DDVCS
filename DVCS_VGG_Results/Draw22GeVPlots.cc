

void Draw22GeVPlots(){
  TCanvas *c1 = new TCanvas("c1", "", 950, 950);
  c1->SetRightMargin(0.04);
  c1->SetLeftMargin(0.12);
  c1->SetLogy();
  
  TLatex lat1;
  lat1.SetNDC();
  lat1.SetTextFont(42);
  lat1.SetTextSize(0.04);

  double Q2 = 2.75;
  double xB = 0.15;
  double tM = -0.3;
  double Qp2 = 0.3;

  
  TGraph gr_DDVCSOnly_tDep_22GeV("DDVCSOnly_22GeV_xsec.dat", "%*lg %lg %lg");
  TGraph gr_DDVCSPlusBH_tDep_22GeV("DDVCSPlusBH_22GeV_xsec.dat", "%*lg %lg %lg");
  TGraph gr_DVCSOnly_tDep_22GeV("DVCSOnly_22GeV_xsec.dat", "%lg %*lg %lg");
  TGraph gr_DVCSPlusBH_tDep_22GeV("DVCSPlusBH_22GeV_xsec.dat", "%lg %*lg %lg");

  TGraph gr_DDVCSOnly_tDep_Set2("DDVCOnly_tDep_Set2_out.dat", "%*lg %lg %lg");
  TGraph gr_DDVCSPlusBH_tDep_Set2("DDVCSPlusBH_tDep_Set2_out.dat", "%*lg %lg %lg");
  TGraph gr_DVCSOnly_tDep_Set2("DVCSOnly_tDep_Set2_out.dat", "%lg %*lg %lg");
  TGraph gr_DVCSPlusBH_tDep_Set2("DVCSPlusBH_tDep_Set2_out.dat", "%lg %*lg %lg");

  gr_DVCSOnly_tDep_Set2.SetMarkerColor(2);
  gr_DVCSOnly_tDep_Set2.SetLineColor(6);
  gr_DVCSOnly_tDep_Set2.SetLineStyle(8);
  gr_DVCSOnly_tDep_Set2.SetLineWidth(3);
  gr_DVCSPlusBH_tDep_Set2.SetMarkerColor(2);
  gr_DVCSPlusBH_tDep_Set2.SetLineColor(6);
  gr_DVCSPlusBH_tDep_Set2.SetLineStyle(1);
  gr_DVCSPlusBH_tDep_Set2.SetLineWidth(3);

  
  gr_DVCSOnly_tDep_22GeV.SetMarkerColor(2);
  gr_DVCSOnly_tDep_22GeV.SetLineColor(2);
  gr_DVCSOnly_tDep_22GeV.SetLineStyle(8);
  gr_DVCSPlusBH_tDep_22GeV.SetMarkerColor(2);
  gr_DVCSPlusBH_tDep_22GeV.SetLineColor(2);
  gr_DVCSPlusBH_tDep_22GeV.SetLineStyle(1);

  TLegend *leg_DVCStDep = new TLegend(0.65, 0.7, 0.95, 0.88);
  leg_DVCStDep->SetBorderSize(0);
  leg_DVCStDep->AddEntry(&gr_DVCSPlusBH_tDep_22GeV, "DVCS + BH: 22 GeV");
  leg_DVCStDep->AddEntry(&gr_DVCSOnly_tDep_22GeV, "DVCS: 22 GeV");
  leg_DVCStDep->AddEntry(&gr_DVCSPlusBH_tDep_Set2, "DVCS + BH: 10.6 GeV");
  leg_DVCStDep->AddEntry(&gr_DVCSOnly_tDep_Set2, "DVCS: 10.6 GeV");

  double phi_LH = 0;
  TMultiGraph *mtgr_DVCS_tDep_2 = new TMultiGraph();
  mtgr_DVCS_tDep_2->Add(&gr_DVCSOnly_tDep_22GeV);
  mtgr_DVCS_tDep_2->Add(&gr_DVCSPlusBH_tDep_22GeV);
  mtgr_DVCS_tDep_2->Add(&gr_DVCSOnly_tDep_Set2);
  mtgr_DVCS_tDep_2->Add(&gr_DVCSPlusBH_tDep_Set2);
  mtgr_DVCS_tDep_2->Draw("Al");
  mtgr_DVCS_tDep_2->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dt dx_{B} d#Phi [nb/GeV^{4}]");
  lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, #Phi=%1.2f deg", Q2, xB, phi_LH));
  leg_DVCStDep->Draw();
  mtgr_DVCS_tDep_2->GetXaxis()->SetLimits(0., 1.5);
  c1->Print("Figs/DVCS_tDep_22GeV.pdf");
  c1->Print("Figs/DVCS_tDep_22GeV.png");
  c1->Print("Figs/DVCS_tDep_22GeV.root");


  gr_DDVCSOnly_tDep_22GeV.SetMarkerColor(4);
  gr_DDVCSOnly_tDep_22GeV.SetLineColor(4);
  gr_DDVCSOnly_tDep_22GeV.SetLineStyle(8);
  gr_DDVCSPlusBH_tDep_22GeV.SetMarkerColor(4);
  gr_DDVCSPlusBH_tDep_22GeV.SetLineColor(4);
  gr_DDVCSPlusBH_tDep_22GeV.SetLineStyle(1);

  gr_DDVCSOnly_tDep_Set2.SetMarkerColor(6);
  gr_DDVCSOnly_tDep_Set2.SetLineColor(6);
  gr_DDVCSOnly_tDep_Set2.SetLineStyle(8);
  gr_DDVCSOnly_tDep_Set2.SetLineWidth(3);
  gr_DDVCSPlusBH_tDep_Set2.SetMarkerColor(6);
  gr_DDVCSPlusBH_tDep_Set2.SetLineColor(6);
  gr_DDVCSPlusBH_tDep_Set2.SetLineStyle(1);
  gr_DDVCSPlusBH_tDep_Set2.SetLineWidth(3);


  TLegend *leg_DDVCStDep = new TLegend(0.65, 0.7, 0.95, 0.88);
  leg_DDVCStDep->SetBorderSize(0);
  leg_DDVCStDep->AddEntry(&gr_DDVCSPlusBH_tDep_22GeV, "DDVCS + BH: 22 GeV");
  leg_DDVCStDep->AddEntry(&gr_DDVCSOnly_tDep_22GeV, "DDVCS: 22 GeV");
  leg_DDVCStDep->AddEntry(&gr_DDVCSPlusBH_tDep_Set2, "DDVCS + BH: 10.6 GeV");
  leg_DDVCStDep->AddEntry(&gr_DDVCSOnly_tDep_Set2, "DDVCS: 10.6 GeV");
  TMultiGraph *mtgr_DDVCS_tDep_2 = new TMultiGraph();
  mtgr_DDVCS_tDep_2->Add(&gr_DDVCSOnly_tDep_22GeV);
  mtgr_DDVCS_tDep_2->Add(&gr_DDVCSPlusBH_tDep_22GeV);
  mtgr_DDVCS_tDep_2->Add(&gr_DDVCSOnly_tDep_Set2);
  mtgr_DDVCS_tDep_2->Add(&gr_DDVCSPlusBH_tDep_Set2);
  mtgr_DDVCS_tDep_2->Draw("Al");
  mtgr_DDVCS_tDep_2->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dt dx_{B} d#Phi dQ^{'2} [nb/GeV^{6}]");
  lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, Q^{'2} = %1.2f GeV^{2}, #Phi=%1.2f deg", Q2, xB, Qp2, phi_LH));
  leg_DDVCStDep->Draw();
  mtgr_DDVCS_tDep_2->GetXaxis()->SetLimits(0., 1.5);
  mtgr_DDVCS_tDep_2->SetMinimum(1.e-8);
  c1->Print("Figs/DDVCS_tDep_22GeV.pdf");
  c1->Print("Figs/DDVCS_tDep_22GeV.png");
  c1->Print("Figs/DDVCS_tDep_22GeV.root");
    
}
