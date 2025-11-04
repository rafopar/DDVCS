#include <cstdlib>
#include <fstream>
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>

using namespace std;



void DrawPlots1(int run) {

    std::map<int, double> m_TotXsec;
    std::map<int, double> m_Eb;
    std::map<int, double> m_totCharge;
    std::map<int, double> m_BeamPol;
    std::map<int, bool> m_isMC;

        std::ifstream inp_config("RunSettings.dat");
    if (!inp_config.is_open()) {
        std::cerr << "Error opening the config file! \n Exiting..." << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(inp_config, line)) {
        size_t commentPos = line.find("#"); // Or line.find("#"); for hash comments

        // If a comment is found, truncate the string at that position
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string entry;
        std::vector<std::string> v_entries;

        while (iss >> entry) {
            v_entries.push_back(entry);
        }

        int arun = std::stoi(v_entries[0]);
        double axSec = std::stof(v_entries[1]);
        double aEb = std::stof(v_entries[2]);
        double aCarge = std::stof(v_entries[3]);
        double aBeamPol = std::stof(v_entries[4]);
        int aisMC = std::stoi(v_entries[5]);

        m_TotXsec[arun] = axSec;
        m_Eb[arun] = aEb;
        m_totCharge[arun] = aCarge;
        m_BeamPol[arun] = aBeamPol;
        m_isMC[arun] = aisMC;
    }



    //bool isMC = data_Set.find("MC") == std::string::npos ? 0 : 1;
    bool isMC = m_isMC[run];

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

    TFile *file_in = new TFile(Form("HistData/Hists_DDVCS_Run_%d.root", run), "Read");

    TF1 *f_GPol3 = new TF1("f_GPol3", " [0] + x*( [1] + x*([2] + x*[3]) ) + [4]*TMath::Gaus(x, [5], [6])", 0.5, 1.5);
    f_GPol3->SetNpx(4500);
    TF1 *f_Pol3 = new TF1("f_Pol3", "[0] + x*( [1] + x*([2] + x*[3]) )", 0.5, 1.5);
    f_Pol3->SetNpx(4500);
    f_Pol3->SetLineColor(6);
    TF1 *f_Gaus = new TF1("f_Gaus", "[0]*TMath::Gaus(x, [1], [2])", 0.5, 1.5);
    f_Gaus->SetNpx(4500);
    f_Gaus->SetLineColor(4);

    //if (strcmp(data_Set.c_str(), "S19") == 0) {
    TH1D *h_Mmis3 = (TH1D*) file_in->Get("h_Mmis3");
    h_Mmis3->SetTitle("; M_{X} [GeV] ");
    h_Mmis3->Draw();
    f_GPol3->SetParameter(4, 100.);
    f_GPol3->SetParameter(5, 1.);
    f_GPol3->SetParameter(6, 0.05);
    h_Mmis3->Fit(f_GPol3, "MeV", "", 0.5, 1.8);
    double pars[7];
    f_GPol3->GetParameters(pars);
    f_Pol3->SetParameters(pars);
    f_Gaus->SetParameters(&pars[4]);
    f_Pol3->Draw("Same");
    f_Gaus->Draw("Same");
    double N_Gaus = f_Gaus->Integral(0.5, 1.5) / h_Mmis3->GetBinWidth(15);
    lat1->DrawLatex(0.2, 0.8, Form("N prot = %1.1f", N_Gaus));
    c1->Print(Form("Figs/MMis3_Run_%d.pdf", run));
    c1->Print(Form("Figs/MMis3_Run_%d.png", run));
    c1->Print(Form("Figs/MMis3_Run_%d.root", run));


    TH2D *h_th_VS_Mmis3 = (TH2D*) file_in->Get(Form("h_th_VS_Mmis3"));
    h_th_VS_Mmis3->SetTitle("; #theta_{Mis} [deg]; M_{Mis} [GeV]");
    h_th_VS_Mmis3->Draw("colz");
    c1->Print(Form("Figs/Mmis_Vs_th3_Run_%d.pdf", run));
    c1->Print(Form("Figs/Mmis_Vs_th3_Run_%d.png", run));
    c1->Print(Form("Figs/Mmis_Vs_th3_Run_%d.root", run));

    const double thMis_Cut = 15; // degree, will use events only above this cut.
    double binThCut = h_th_VS_Mmis3->GetXaxis()->FindBin(thMis_Cut);
    TH1D *h_Mmis3_ThMisCut = (TH1D*) h_th_VS_Mmis3->ProjectionY("h_Mmis3_ThMisCut", binThCut, h_th_VS_Mmis3->GetNbinsX());
    h_Mmis3_ThMisCut->Draw();
    h_Mmis3_ThMisCut->Fit(f_GPol3, "MeV", "", 0.5, 1.8);
    f_GPol3->GetParameters(pars);
    f_Pol3->SetParameters(pars);
    f_Gaus->SetParameters(&pars[4]);
    f_Pol3->Draw("Same");
    f_Gaus->Draw("Same");
    N_Gaus = f_Gaus->Integral(0.5, 1.5) / h_Mmis3_ThMisCut->GetBinWidth(15);
    lat1->DrawLatex(0.2, 0.8, Form("N prot = %1.1f", N_Gaus));
    c1->Print(Form("Figs/Mmis3_ThMisCut_Run_%d.pdf", run));
    c1->Print(Form("Figs/Mmis3_ThMisCut_Run_%d.png", run));
    c1->Print(Form("Figs/Mmis3_ThMisCut_Run_%d.root", run));

    TH2D *h_th_VS_Mmis4 = (TH2D*) file_in->Get(Form("h_th_VS_Mmis4"));
    h_th_VS_Mmis4->SetTitle("; #theta_{Mis} [deg]; M_{Mis} [GeV]");
    h_th_VS_Mmis4->Draw("colz");
    c1->Print(Form("Figs/Mmis_Vs_th4_Run_%d.pdf", run));
    c1->Print(Form("Figs/Mmis_Vs_th4_Run_%d.png", run));
    c1->Print(Form("Figs/Mmis_Vs_th4_Run_%d.root", run));

    binThCut = h_th_VS_Mmis4->GetXaxis()->FindBin(thMis_Cut);
    TH1D *h_Mmis4_ThMisCut = (TH1D*) h_th_VS_Mmis4->ProjectionY("h_Mmis4_ThMisCut", binThCut, h_th_VS_Mmis4->GetNbinsX());
    h_Mmis4_ThMisCut->Draw();
    h_Mmis4_ThMisCut->Fit(f_GPol3, "MeV", "", 0.5, 1.8);
    f_GPol3->GetParameters(pars);
    f_Pol3->SetParameters(pars);
    f_Gaus->SetParameters(&pars[4]);
    f_Pol3->Draw("Same");
    f_Gaus->Draw("Same");
    N_Gaus = f_Gaus->Integral(0.5, 1.5) / h_Mmis4_ThMisCut->GetBinWidth(15);
    lat1->DrawLatex(0.2, 0.8, Form("N prot = %1.1f", N_Gaus));
    c1->Print(Form("Figs/Mmis4_ThMisCut_Run_%d.pdf", run));
    c1->Print(Form("Figs/Mmis4_ThMisCut_Run_%d.png", run));
    c1->Print(Form("Figs/Mmis4_ThMisCut_Run_%d.root", run));

    TH1D *h_Mmis1 = (TH1D*) file_in->Get("h_Mmis1");
    h_Mmis1->SetTitle("; M_{X} [GeV] ");
    h_Mmis1->Draw("hist");
    cout << "isMC = " << isMC << endl;
    if (isMC) {
        lat1->DrawLatex(0.55, 0.6, Form("Tot events = %1.1f", h_Mmis1->Integral()));
    }
    c1->Print(Form("Figs/MMis1_Run_%d.pdf", run));
    c1->Print(Form("Figs/MMis1_Run_%d.png", run));
    c1->Print(Form("Figs/MMis1_Run_%d.root", run));

    TH2D *h_Mmis_PMis1 = (TH2D*)file_in->Get("h_Mmis_PMis1");
    h_Mmis_PMis1->SetTitle("; P_{Mis} [GeV]; M_{Mis} [GeV]");
    h_Mmis_PMis1->Draw("colz");
    c1->Print(Form("Figs/Mmis_PMis1_Run_%d.pdf", run));
    c1->Print(Form("Figs/Mmis_PMis1_Run_%d.png", run));
    c1->Print(Form("Figs/Mmis_PMis1_Run_%d.root", run));
    
    TH2D *h_Mmis_PMis3 = (TH2D*)file_in->Get("h_Mmis_PMis3");
    h_Mmis_PMis3->SetTitle("; P_{Mis} [GeV]; M_{Mis} [GeV]");
    h_Mmis_PMis3->Draw("colz");
    c1->Print(Form("Figs/Mmis_PMis3_Run_%d.pdf", run));
    c1->Print(Form("Figs/Mmis_PMis3_Run_%d.png", run));
    c1->Print(Form("Figs/Mmis_PMis3_Run_%d.root", run));
    
    
    if (isMC) {
        TH1D *h_Mmis3 = (TH1D*) file_in->Get("h_Mmis3");
        h_Mmis3->SetTitle("; M_{X} [GeV] ");
        h_Mmis3->Draw("hist");
        cout << "isMC = " << isMC << endl;
        lat1->DrawLatex(0.35, 0.6, Form("Tot events = %1.1f", h_Mmis3->Integral()));
        double peak_right_edge = 1.1;
        int bin_peak_edge = h_Mmis3->FindBin(peak_right_edge); // after 1.1 GeV the peak finished and all is essentially the radiative tail.
        double ev_under_peak = h_Mmis3->Integral(1, bin_peak_edge);
        lat1->DrawLatex(0.35, 0.5, Form("Events belw %1.2f GeV is %1.1f", peak_right_edge, ev_under_peak));
        c1->Print(Form("Figs/MMis3_Run_%d.pdf", run));
        c1->Print(Form("Figs/MMis3_Run_%d.png", run));
        c1->Print(Form("Figs/MMis3_Run_%d.root", run));
    }

    TH1D *h_vt_Diff_em1 = (TH1D*) file_in->Get(Form("h_vt_Diff_em1"));
    h_vt_Diff_em1->SetTitle("; t_e^{-}_{1} - t_e^{-}_{2} [ns]");
    h_vt_Diff_em1->Draw("hist");
    c1->Print(Form("Figs/vt_diff_em1_Run_%d.pdf", run));
    c1->Print(Form("Figs/vt_diff_em1_Run_%d.png", run));
    c1->Print(Form("Figs/vt_diff_em1_Run_%d.root", run));

    TH2D *h_vt_Diff_emep1 = (TH2D*) file_in->Get("h_vt_Diff_emep1");
    h_vt_Diff_emep1->SetTitle("; t_e^{-}_{1} - t_e^{-}_{2} [ns]; t_e^{-}_{1} - t_e^{+} [ns]");
    h_vt_Diff_emep1->Draw("colz");
    c1->Print(Form("Figs/vt_diff_em1_ep_Run_%d.pdf", run));
    c1->Print(Form("Figs/vt_diff_em1_ep_Run_%d.png", run));
    c1->Print(Form("Figs/vt_diff_em1_ep_Run_%d.root", run));

    TH1D *h_vz_em1 = (TH1D*) file_in->Get(Form("h_vz_em1"));
    h_vz_em1->SetTitle("; vz e^{-} [mm]");
    h_vz_em1->Draw("hist");
    c1->Print(Form("Figs/vz_em1_Run_%d.pdf", run));
    c1->Print(Form("Figs/vz_em1_Run_%d.png", run));
    c1->Print(Form("Figs/vz_em1_Run_%d.root", run));

    TH1D *h_vz_ep1 = (TH1D*) file_in->Get(Form("h_vz_ep1"));
    h_vz_ep1->SetTitle("; vz e^{+} [mm]");
    h_vz_ep1->Draw("hist");
    c1->Print(Form("Figs/vz_ep1_Run_%d.pdf", run));
    c1->Print(Form("Figs/vz_ep1_Run_%d.png", run));
    c1->Print(Form("Figs/vz_ep1_Run_%d.root", run));

    TH2D *h_vz_Diff_emep1 = (TH2D*) file_in->Get("h_vz_Diff_emep1");
    h_vz_Diff_emep1->SetTitle("; vz_e^{-}_{1} - vz_e^{-}_{2} [cm]; vz_e^{-}_{1} - vz_e^{+} [cm]");
    h_vz_Diff_emep1->Draw("colz");
    c1->Print(Form("Figs/vz_diff_em1_ep_Run_%d.pdf", run));
    c1->Print(Form("Figs/vz_diff_em1_ep_Run_%d.png", run));
    c1->Print(Form("Figs/vz_diff_em1_ep_Run_%d.root", run));

    TH2D *h_Minv12_1 = (TH2D*) file_in->Get("h_Minv12_1");
    h_Minv12_1->SetTitle(";M(e^{-}_{1}e^{+}) [GeV]; M(e^{-}_{2}e^{+}) [GeV]");
    h_Minv12_1->SetStats(0);
    h_Minv12_1->Draw("colz");
    line1->DrawLine(M_eeCut, M_eeCut, M_eeCut, 2.5);
    line1->DrawLine(M_eeCut, M_eeCut, 2.5, M_eeCut);
    c1->Print(Form("Figs/Minv1_Minv2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Minv1_Minv2_Run_%d.png", run));
    c1->Print(Form("Figs/Minv1_Minv2_Run_%d.root", run));

    TH2D *h_Minv12_2 = (TH2D*) file_in->Get("h_Minv12_2");
    h_Minv12_2->SetTitle(";M(e^{-}_{1}e^{+}) [GeV]; M(e^{-}_{2}e^{+}) [GeV]");
    h_Minv12_2->SetStats(0);
    h_Minv12_2->Draw("colz");
    line1->DrawLine(M_eeCut, M_eeCut, M_eeCut, 2.5);
    line1->DrawLine(M_eeCut, M_eeCut, 2.5, M_eeCut);
    c1->Print(Form("Figs/Minv1_Minv2_2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Minv1_Minv2_2_Run_%d.png", run));
    c1->Print(Form("Figs/Minv1_Minv2_2_Run_%d.root", run));

    c1->SetLogy();
    TH1D *h_vt_Diff_em2 = (TH1D*) file_in->Get("h_vt_Diff_em2");
    h_vt_Diff_em2->SetTitle("; t_e^{-}_{1} - t_e^{-}_{2} [ns]");
    h_vt_Diff_em2->Draw("hist");
    line1->DrawLine(dtCut, 0., dtCut, h_vt_Diff_em2->GetMaximum());
    line1->DrawLine(-dtCut, 0., -dtCut, h_vt_Diff_em2->GetMaximum());
    c1->Print(Form("Figs/vt_Diff_em2_Run_%d.pdf", run));
    c1->Print(Form("Figs/vt_Diff_em2_Run_%d.png", run));
    c1->Print(Form("Figs/vt_Diff_em2_Run_%d.root", run));

    c1->SetLogy(0);
    
    TH2D *h_th_P_em1 = (TH2D*)file_in->Get("h_th_P_em1");
    h_th_P_em1->SetTitle("; P_{e^{-}} [GeV]; #theta_{ e^{-} } [deg]");
    h_th_P_em1->Draw("colz");
    c1->Print(Form("Figs/th_P_em1_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em1_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em1_Run_%d.root", run));

    TH2D *h_th_P_ep1 = (TH2D*)file_in->Get("h_th_P_ep1");
    h_th_P_ep1->SetTitle("; P_{e^{+}} [GeV]; #theta_{ e^{+} } [deg]");
    h_th_P_ep1->Draw("colz");
    c1->Print(Form("Figs/th_P_ep1_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_ep1_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_ep1_Run_%d.root", run));

    TH2D *h_th_P_em2 = (TH2D*)file_in->Get("h_th_P_em2");
    h_th_P_em2->SetTitle("; P_{e^{-}} [GeV]; #theta_{ e^{-} } [deg]");
    h_th_P_em2->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em2_Run_%d.root", run));
    
    TH2D *h_th_P_em1_2 = (TH2D*)file_in->Get("h_th_P_em1_2");
    h_th_P_em1_2->SetTitle("; P_{e^{-}_{1}} [GeV]; #theta_{ e^{-}_{1} } [deg]");
    h_th_P_em1_2->Draw("colz");
    c1->Print(Form("Figs/th_P_em1_2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em1_2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em1_2_Run_%d.root", run));

    TH2D *h_th_P_em2_2 = (TH2D*)file_in->Get("h_th_P_em2_2");
    h_th_P_em2_2->SetTitle("; P_{e^{-}_{2}} [GeV]; #theta_{ e^{-}_{2} } [deg]");
    h_th_P_em2_2->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em2_2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em2_2_Run_%d.root", run));

    TH2D *h_th_P_ep2 = (TH2D*)file_in->Get("h_th_P_ep2");
    h_th_P_ep2->SetTitle("; P_{e^{+}} [GeV]; #theta_{ e^{+} } [deg]");
    h_th_P_ep2->Draw("colz");
    c1->Print(Form("Figs/th_P_ep2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_ep2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_ep2_Run_%d.root", run));

    TH2D *h_th_P_em3 = (TH2D*)file_in->Get("h_th_P_em3");
    h_th_P_em3->SetTitle("; P_{e^{-}} [GeV]; #theta_{ e^{-} } [deg]");
    h_th_P_em3->Draw("colz");
    c1->Print(Form("Figs/th_P_em3_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em3_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em3_Run_%d.root", run));
    
    TH2D *h_th_P_em1_3 = (TH2D*)file_in->Get("h_th_P_em1_3");
    h_th_P_em1_3->SetTitle("; P_{e^{-}_{1}} [GeV]; #theta_{ e^{-}_{1} } [deg]");
    h_th_P_em1_3->Draw("colz");
    c1->Print(Form("Figs/th_P_em1_3_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em1_3_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em1_3_Run_%d.root", run));

    TH2D *h_th_P_em2_3 = (TH2D*)file_in->Get("h_th_P_em2_3");
    h_th_P_em2_3->SetTitle("; P_{e^{-}_{2}} [GeV]; #theta_{ e^{-}_{2} } [deg]");
    h_th_P_em2_3->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_3_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em2_3_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em2_3_Run_%d.root", run));

    TH2D *h_th_P_ep3 = (TH2D*)file_in->Get("h_th_P_ep3");
    h_th_P_ep3->SetTitle("; P_{e^{+}} [GeV]; #theta_{ e^{+} } [deg]");
    h_th_P_ep3->Draw("colz");
    c1->Print(Form("Figs/th_P_ep3_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_ep3_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_ep3_Run_%d.root", run));

    
    TH2D *h_th_P_em4 = (TH2D*)file_in->Get("h_th_P_em4");
    h_th_P_em4->SetTitle("; P_{e^{-}} [GeV]; #theta_{ e^{-} } [deg]");
    h_th_P_em4->Draw("colz");
    c1->Print(Form("Figs/th_P_em4_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em4_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em4_Run_%d.root", run));
    
    TH2D *h_th_P_em1_4 = (TH2D*)file_in->Get("h_th_P_em1_4");
    h_th_P_em1_4->SetTitle("; P_{e^{-}_{1}} [GeV]; #theta_{ e^{-}_{1} } [deg]");
    h_th_P_em1_4->Draw("colz");
    c1->Print(Form("Figs/th_P_em1_4_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em1_4_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em1_4_Run_%d.root", run));

    TH2D *h_th_P_em2_4 = (TH2D*)file_in->Get("h_th_P_em2_4");
    h_th_P_em2_4->SetTitle("; P_{e^{-}_{2}} [GeV]; #theta_{ e^{-}_{2} } [deg]");
    h_th_P_em2_4->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_4_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em2_4_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em2_4_Run_%d.root", run));

    TH2D *h_th_P_ep4 = (TH2D*)file_in->Get("h_th_P_ep4");
    h_th_P_ep4->SetTitle("; P_{e^{+}} [GeV]; #theta_{ e^{+} } [deg]");
    h_th_P_ep4->Draw("colz");
    c1->Print(Form("Figs/th_P_ep4_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_ep4_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_ep4_Run_%d.root", run));

    
    TH2D *h_th_P_em5 = (TH2D*)file_in->Get("h_th_P_em5");
    h_th_P_em5->SetTitle("; P_{e^{-}} [GeV]; #theta_{ e^{-} } [deg]");
    h_th_P_em5->Draw("colz");
    c1->Print(Form("Figs/th_P_em5_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em5_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em5_Run_%d.root", run));
    
    TH2D *h_th_P_em1_5 = (TH2D*)file_in->Get("h_th_P_em1_5");
    h_th_P_em1_5->SetTitle("; P_{e^{-}_{1}} [GeV]; #theta_{ e^{-}_{1} } [deg]");
    h_th_P_em1_5->Draw("colz");
    c1->Print(Form("Figs/th_P_em1_5_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em1_5_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em1_5_Run_%d.root", run));

    TH2D *h_th_P_em2_5 = (TH2D*)file_in->Get("h_th_P_em2_5");
    h_th_P_em2_5->SetTitle("; P_{e^{-}_{2}} [GeV]; #theta_{ e^{-}_{2} } [deg]");
    h_th_P_em2_5->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_5_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em2_5_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em2_5_Run_%d.root", run));

    TH2D *h_th_P_ep5 = (TH2D*)file_in->Get("h_th_P_ep5");
    h_th_P_ep5->SetTitle("; P_{e^{+}} [GeV]; #theta_{ e^{+} } [deg]");
    h_th_P_ep5->Draw("colz");
    c1->Print(Form("Figs/th_P_ep5_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_ep5_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_ep5_Run_%d.root", run));
    
    TH1D *h_th_em5 = (TH1D*)h_th_P_em5->ProjectionY("h_th_em5", 1, h_th_P_em5->GetNbinsX() );
    h_th_em5->Draw();
    c1->Print(Form("Figs/th_em5_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_em5_Run_%d.png", run));
    c1->Print(Form("Figs/th_em5_Run_%d.root", run));

    TH1D *h_P_em5 = (TH1D*)h_th_P_em5->ProjectionX("h_P_em5", 1, h_th_P_em5->GetNbinsY() );
    h_P_em5->Draw();
    c1->Print(Form("Figs/P_em5_Run_%d.pdf", run));
    c1->Print(Form("Figs/P_em5_Run_%d.png", run));
    c1->Print(Form("Figs/P_em5_Run_%d.root", run));
  
    TH1D *h_th_ep5 = (TH1D*)h_th_P_ep5->ProjectionY("h_th_ep5", 1, h_th_P_ep5->GetNbinsX() );
    h_th_ep5->Draw();
    c1->Print(Form("Figs/th_ep5_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_ep5_Run_%d.png", run));
    c1->Print(Form("Figs/th_ep5_Run_%d.root", run));

    TH1D *h_P_ep5 = (TH1D*)h_th_P_ep5->ProjectionX("h_P_ep5", 1, h_th_P_ep5->GetNbinsY() );
    h_P_ep5->Draw();
    c1->Print(Form("Figs/P_ep5_Run_%d.pdf", run));
    c1->Print(Form("Figs/P_ep5_Run_%d.png", run));
    c1->Print(Form("Figs/P_ep5_Run_%d.root", run));

    const double PCal_VMinCut = 10;
    const double PCal_WMinCut = 10;
    TH1D *h_lVPCal_em1 = (TH1D*)file_in->Get("h_lVPCal_em1");
    h_lVPCal_em1->SetTitle("; PCal V (e^{-}) [cm]");
    h_lVPCal_em1->Draw("hist");
    line1->DrawLine(PCal_VMinCut, 0, PCal_VMinCut, h_lVPCal_em1->GetMaximum() );
    c1->Print(Form("Figs/lVPCal_em1_Run_%d.pdf", run));
    c1->Print(Form("Figs/lVPCal_em1_Run_%d.png", run));
    c1->Print(Form("Figs/lVPCal_em1_Run_%d.root", run));

    TH1D *h_lWPCal_em1 = (TH1D*)file_in->Get("h_lWPCal_em1");
    h_lWPCal_em1->SetTitle("; PCal W (e^{-}) [cm]");
    h_lWPCal_em1->Draw("hist");
    line1->DrawLine(PCal_WMinCut, 0, PCal_WMinCut, h_lWPCal_em1->GetMaximum() );
    c1->Print(Form("Figs/lWPCal_em1_Run_%d.pdf", run));
    c1->Print(Form("Figs/lWPCal_em1_Run_%d.png", run));
    c1->Print(Form("Figs/lWPCal_em1_Run_%d.root", run));

    TH1D *h_lVPCal_ep1 = (TH1D*)file_in->Get("h_lVPCal_ep1");
    h_lVPCal_ep1->SetTitle("; PCal V (e^{+}) [cm]");
    h_lVPCal_ep1->Draw("hist");
    line1->DrawLine(PCal_VMinCut, 0, PCal_VMinCut, h_lVPCal_ep1->GetMaximum() );
    c1->Print(Form("Figs/lVPCal_ep1_Run_%d.pdf", run));
    c1->Print(Form("Figs/lVPCal_ep1_Run_%d.png", run));
    c1->Print(Form("Figs/lVPCal_ep1_Run_%d.root", run));

    TH1D *h_lWPCal_ep1 = (TH1D*)file_in->Get("h_lWPCal_ep1");
    h_lWPCal_ep1->SetTitle("; PCal W (e^{+}) [cm]");
    h_lWPCal_ep1->Draw("hist");
    line1->DrawLine(PCal_WMinCut, 0, PCal_WMinCut, h_lWPCal_ep1->GetMaximum() );
    c1->Print(Form("Figs/lWPCal_ep1_Run_%d.pdf", run));
    c1->Print(Form("Figs/lWPCal_ep1_Run_%d.png", run));
    c1->Print(Form("Figs/lWPCal_ep1_Run_%d.root", run));


    c1->SetLogz();
    if (isMC) {
        TH2D *h_MC_Memep_12_1 = (TH2D*) file_in->Get("h_MC_Memep_12_1");
        h_MC_Memep_12_1->SetTitle(";M(e^{-}_{1}e^{+}) [GeV]; M(e^{-}_{2}e^{+}) [GeV]");
        h_MC_Memep_12_1->Draw("colz");
        c1->Print(Form("Figs/MC_Minv1_Minv2_Run_%d.pdf", run));
        c1->Print(Form("Figs/MC_Minv1_Minv2_Run_%d.png", run));
        c1->Print(Form("Figs/MC_Minv1_Minv2_Run_%d.root", run));
    }
}