

void DrawPlots1() {

    //std::string data_type = "S19";
    //std::string data_type = "GRAPE_Run5";
    std::string data_type = "GRAPE_Run6";
    //std::string data_type = "GRAPE_Run8";

    bool isMC = strcmp(data_type.c_str(), "S19") == 0 ? 0 : 1;

    const double M_eeCut = 0.25;
    const double dtCut = 2.5;
    
    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetRightMargin(0.04);
    c1->SetTopMargin(0.03);
    c1->SetLeftMargin(0.14);
    
    TLine *line1 = new TLine();
    line1->SetLineColor(2);
    line1->SetLineWidth(2);

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);

    TFile *file_in = new TFile(Form("Hists_DDVCS_%s.root", data_type.c_str()), "Read");

    TF1 *f_GPol3 = new TF1("f_GPol3", " [0] + x*( [1] + x*([2] + x*[3]) ) + [4]*TMath::Gaus(x, [5], [6])", 0.5, 1.5);
    f_GPol3->SetNpx(4500);
    TF1 *f_Pol3 = new TF1("f_Pol3", "[0] + x*( [1] + x*([2] + x*[3]) )", 0.5, 1.5);
    f_Pol3->SetNpx(4500);
    TF1 *f_Gaus = new TF1("f_Gaus", "[0]*TMath::Gaus(x, [1], [2])", 0.5, 1.5);
    f_Gaus->SetNpx(4500);
    f_Gaus->SetLineColor(4);

    if (strcmp(data_type.c_str(), "S19") == 0) {
        TH1D *h_Mmis3 = (TH1D*) file_in->Get("h_Mmis3");
        h_Mmis3->SetTitle("; M_{X} [GeV] ");
        h_Mmis3->Draw();
        f_GPol3->SetParameter(4, 20.);
        f_GPol3->SetParameter(5, 1.);
        f_GPol3->SetParameter(6, 0.1);
        h_Mmis3->Fit(f_GPol3, "MeV", "", 0.5, 1.8);
        double pars[7];
        f_GPol3->GetParameters(pars);
        f_Pol3->SetParameters(pars);
        f_Gaus->SetParameters(&pars[4]);
        f_Pol3->Draw("Same");
        f_Gaus->Draw("Same");
        double N_Gaus = f_Gaus->Integral(0.5, 1.5) / h_Mmis3->GetBinWidth(15);
        lat1->DrawLatex(0.2, 0.8, Form("N prot = %1.1f", N_Gaus));
        c1->Print(Form("Figs/MMis3_%s.pdf", data_type.c_str()));
        c1->Print(Form("Figs/MMis3_%s.png", data_type.c_str()));
        c1->Print(Form("Figs/MMis3_%s.root", data_type.c_str()));
    }

    TH1D *h_Mmis1 = (TH1D*) file_in->Get("h_Mmis1");
    h_Mmis1->SetTitle("; M_{X} [GeV] ");
    h_Mmis1->Draw("hist");
    cout<<"isMC = "<<isMC<<endl;
    if (isMC) {       
        lat1->DrawLatex(0.55, 0.6, Form("Tot events = %1.1f", h_Mmis1->Integral() ));
    }
    c1->Print(Form("Figs/MMis1_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/MMis1_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/MMis1_%s.root", data_type.c_str()));

    TH1D *h_Mmis3 = (TH1D*) file_in->Get("h_Mmis3");
    h_Mmis3->SetTitle("; M_{X} [GeV] ");
    h_Mmis3->Draw("hist");
    cout<<"isMC = "<<isMC<<endl;
    if (isMC) {       
        lat1->DrawLatex(0.35, 0.6, Form("Tot events = %1.1f", h_Mmis3->Integral() ));
        double peak_right_edge = 1.1;
        int bin_peak_edge = h_Mmis3->FindBin(peak_right_edge);  // after 1.1 GeV the peak finished and all is essentially the radiative tail.
        double ev_under_peak = h_Mmis3->Integral(1, bin_peak_edge);
        lat1->DrawLatex(0.35, 0.5, Form("Events belw %1.2f GeV is %1.1f", peak_right_edge, ev_under_peak));
    }
    c1->Print(Form("Figs/MMis3_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/MMis3_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/MMis3_%s.root", data_type.c_str()));

    TH1D *h_vt_Diff_em1 = (TH1D*)file_in->Get(Form("h_vt_Diff_em1"));
    h_vt_Diff_em1->SetTitle("; t_e^{-}_{1} - t_e^{-}_{2} [ns]");
    h_vt_Diff_em1->Draw("hist");
    c1->Print(Form("Figs/vt_diff_em1_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/vt_diff_em1_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/vt_diff_em1_%s.root", data_type.c_str()));
    
    TH2D *h_vt_Diff_emep1 = (TH2D*)file_in->Get("h_vt_Diff_emep1");
    h_vt_Diff_emep1->SetTitle("; t_e^{-}_{1} - t_e^{-}_{2} [ns]; t_e^{-}_{1} - t_e^{+} [ns]");
    h_vt_Diff_emep1->Draw("colz");
    c1->Print(Form("Figs/vt_diff_em1_ep_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/vt_diff_em1_ep_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/vt_diff_em1_ep_%s.root", data_type.c_str()));

    TH2D *h_vz_Diff_emep1 = (TH2D*)file_in->Get("h_vz_Diff_emep1");
    h_vz_Diff_emep1->SetTitle("; vz_e^{-}_{1} - vz_e^{-}_{2} [cm]; vz_e^{-}_{1} - vz_e^{+} [cm]");
    h_vz_Diff_emep1->Draw("colz");
    c1->Print(Form("Figs/vz_diff_em1_ep_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/vz_diff_em1_ep_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/vz_diff_em1_ep_%s.root", data_type.c_str()));
    
    TH2D *h_Minv12_1 = (TH2D*)file_in->Get("h_Minv12_1");
    h_Minv12_1->SetTitle(";M(e^{-}_{1}e^{+}) [GeV]; M(e^{-}_{2}e^{+}) [GeV]");
    h_Minv12_1->SetStats(0);
    h_Minv12_1->Draw("colz");
    line1->DrawLine( M_eeCut, M_eeCut, M_eeCut, 2.5 );
    line1->DrawLine( M_eeCut, M_eeCut, 2.5, M_eeCut );
    c1->Print(Form("Figs/Minv1_Minv2_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/Minv1_Minv2_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/Minv1_Minv2_%s.root", data_type.c_str()));
    
    c1->SetLogy();
    TH1D *h_vt_Diff_em2 = (TH1D*)file_in->Get("h_vt_Diff_em2");
    h_vt_Diff_em2->SetTitle("; t_e^{-}_{1} - t_e^{-}_{2} [ns]");
    h_vt_Diff_em2->Draw();
    line1->DrawLine(dtCut, 0., dtCut, h_vt_Diff_em2->GetMaximum());
    line1->DrawLine(-dtCut, 0., -dtCut, h_vt_Diff_em2->GetMaximum());
    c1->Print(Form("Figs/vt_Diff_em2_%s.pdf", data_type.c_str()));
    c1->Print(Form("Figs/vt_Diff_em2_%s.png", data_type.c_str()));
    c1->Print(Form("Figs/vt_Diff_em2_%s.root", data_type.c_str()));
    
 
    c1->SetLogy(0);
    c1->SetLogz();
    if( isMC ){
        TH2D *h_MC_Memep_12_1 = (TH2D*)file_in->Get("h_MC_Memep_12_1");
        h_MC_Memep_12_1->SetTitle(";M(e^{-}_{1}e^{+}) [GeV]; M(e^{-}_{2}e^{+}) [GeV]");
        h_MC_Memep_12_1->Draw("colz");
        c1->Print(Form("Figs/MC_Minv1_Minv2_%s.pdf", data_type.c_str() ));
        c1->Print(Form("Figs/MC_Minv1_Minv2_%s.png", data_type.c_str() ));
        c1->Print(Form("Figs/MC_Minv1_Minv2_%s.root", data_type.c_str() ));
    }
}