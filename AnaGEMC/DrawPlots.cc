/* 
 * File:   DrawPlots.cc
 * Author: rafopar
 *
 * Created on October 2, 2024, 8:13 PM
 */

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>

using namespace std;

template <typename T> void FormatHist(T h);

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Please provide the run number" << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);
    lat1->SetTextColor(4);
    lat1->SetTextSize(0.04);

    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    int run = atoi(argv[1]);

    int m_cols[2][2];
    m_cols[0][0] = 4;
    m_cols[0][1] = 2;
    m_cols[1][0] = 6;
    m_cols[1][1] = 3;
    
    TCanvas c1("c1", "", 950, 950);
    c1.SetTopMargin(0.04);
    c1.SetRightMargin(0.04);
    c1.SetLeftMargin(0.13);
    c1.SetBottomMargin(0.13);

    TFile file_in(Form("AnaDDVCS_Run_%d.root", run), "Read");

    TH2D *h_MC_th_P_mup_det1 = (TH2D*) file_in.Get("h_MC_th_P_mup_det1");
    h_MC_th_P_mup_det1->SetTitle("; P(#mu^{+}) [GeV]; #theta(#mu^{+}) [deg]");
    FormatHist(h_MC_th_P_mup_det1);
    h_MC_th_P_mup_det1->Draw();
    lat1->DrawLatex(0.55, 0.9, Form("Tot = %3.2e", h_MC_th_P_mup_det1->Integral()));
    c1.Print(Form("Figs/MC_th_P_mup_det1_Run_%d.pdf", run));
    c1.Print(Form("Figs/MC_th_P_mup_det1_Run_%d.png", run));
    c1.Print(Form("Figs/MC_th_P_mup_det1_Run_%d.root", run));

    TH2D *h_MC_th_P_mup_miss1 = (TH2D*) file_in.Get("h_MC_th_P_mup_miss1");
    h_MC_th_P_mup_miss1->SetTitle("; P(#mu^{+}) [GeV]; #theta(#mu^{+}) [deg]");
    FormatHist(h_MC_th_P_mup_miss1);
    h_MC_th_P_mup_miss1->Draw();
    lat1->DrawLatex(0.55, 0.9, Form("Tot = %3.2e", h_MC_th_P_mup_miss1->Integral()));
    c1.Print(Form("Figs/MC_th_P_mup_miss1_Run_%d.pdf", run));
    c1.Print(Form("Figs/MC_th_P_mup_miss1_Run_%d.png", run));
    c1.Print(Form("Figs/MC_th_P_mup_miss1_Run_%d.root", run));

    TH2D *h_MC_th_P_mum_det1 = (TH2D*) file_in.Get("h_MC_th_P_mum_det1");
    h_MC_th_P_mum_det1->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) [deg]");
    FormatHist(h_MC_th_P_mum_det1);
    h_MC_th_P_mum_det1->Draw();
    lat1->DrawLatex(0.55, 0.9, Form("Tot = %3.2e", h_MC_th_P_mum_det1->Integral()));
    c1.Print(Form("Figs/MC_th_P_mum_det1_Run_%d.pdf", run));
    c1.Print(Form("Figs/MC_th_P_mum_det1_Run_%d.png", run));
    c1.Print(Form("Figs/MC_th_P_mum_det1_Run_%d.root", run));

    TH2D *h_MC_th_P_mum_miss1 = (TH2D*) file_in.Get("h_MC_th_P_mum_miss1");
    h_MC_th_P_mum_miss1->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) [deg]");
    FormatHist(h_MC_th_P_mum_miss1);
    h_MC_th_P_mum_miss1->Draw();
    lat1->DrawLatex(0.55, 0.9, Form("Tot = %3.2e", h_MC_th_P_mum_miss1->Integral()));
    c1.Print(Form("Figs/MC_th_P_mum_miss1_Run_%d.pdf", run));
    c1.Print(Form("Figs/MC_th_P_mum_miss1_Run_%d.png", run));
    c1.Print(Form("Figs/MC_th_P_mum_miss1_Run_%d.root", run));

    TH2D *h_MC_th_P_em_det1 = (TH2D*) file_in.Get("h_MC_th_P_em_det1");
    h_MC_th_P_em_det1->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-}) [deg]");
    FormatHist(h_MC_th_P_em_det1);
    h_MC_th_P_em_det1->Draw();
    lat1->DrawLatex(0.55, 0.9, Form("Tot = %3.2e", h_MC_th_P_em_det1->Integral()));
    c1.Print(Form("Figs/MC_th_P_em_det1_Run_%d.pdf", run));
    c1.Print(Form("Figs/MC_th_P_em_det1_Run_%d.png", run));
    c1.Print(Form("Figs/MC_th_P_em_det1_Run_%d.root", run));

    TH2D *h_MC_th_P_em_miss1 = (TH2D*) file_in.Get("h_MC_th_P_em_miss1");
    h_MC_th_P_em_miss1->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-}) [deg]");
    FormatHist(h_MC_th_P_em_miss1);
    h_MC_th_P_em_miss1->Draw();
    lat1->DrawLatex(0.55, 0.9, Form("Tot = %3.2e", h_MC_th_P_em_miss1->Integral()));
    c1.Print(Form("Figs/MC_th_P_em_miss1_Run_%d.pdf", run));
    c1.Print(Form("Figs/MC_th_P_em_miss1_Run_%d.png", run));
    c1.Print(Form("Figs/MC_th_P_em_miss1_Run_%d.root", run));

    TH2D *h_th_P_mup1 = (TH2D*) file_in.Get("h_th_P_mup1");
    h_th_P_mup1->SetTitle("; P_{#mu^{+}} [GeV]; #theta_{#mu^{+}} ");
    FormatHist(h_th_P_mup1);
    h_th_P_mup1->Draw();
    c1.Print(Form("Figs/th_P_mup1_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mup1_%d.png", run));
    c1.Print(Form("Figs/th_P_mup1_%d.root", run));

    TH2D *h_th_P_mup2 = (TH2D*) file_in.Get("h_th_P_mup2");
    h_th_P_mup2->SetTitle("; P_{#mu^{+}} [GeV]; #theta_{#mu^{+}} ");
    FormatHist(h_th_P_mup2);
    h_th_P_mup2->Draw();
    c1.Print(Form("Figs/th_P_mup2_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mup2_%d.png", run));
    c1.Print(Form("Figs/th_P_mup2_%d.root", run));

    TH2D *h_th_P_mum1 = (TH2D*) file_in.Get("h_th_P_mum1");
    h_th_P_mum1->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) ");
    FormatHist(h_th_P_mum1);
    h_th_P_mum1->Draw();
    c1.Print(Form("Figs/th_P_mum1_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mum1_%d.png", run));
    c1.Print(Form("Figs/th_P_mum1_%d.root", run));

    TH2D *h_th_P_mum2 = (TH2D*) file_in.Get("h_th_P_mum2");
    h_th_P_mum2->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) ");
    FormatHist(h_th_P_mum2);
    h_th_P_mum2->Draw();
    c1.Print(Form("Figs/th_P_mum2_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mum2_%d.png", run));
    c1.Print(Form("Figs/th_P_mum2_%d.root", run));

    TH2D *h_MC_th_P_em1 = (TH2D*) file_in.Get("h_MC_th_P_em1");
    h_MC_th_P_em1->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-}) ");
    FormatHist(h_MC_th_P_em1);
    h_MC_th_P_em1->Draw();
    c1.Print(Form("Figs/th_P_em1_%d.pdf", run));
    c1.Print(Form("Figs/th_P_em1_%d.png", run));
    c1.Print(Form("Figs/th_P_em1_%d.root", run));

    TH2D *h_MC_th_P_em2 = (TH2D*) file_in.Get("h_MC_th_P_em2");
    h_MC_th_P_em2->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-}) ");
    FormatHist(h_MC_th_P_em2);
    h_MC_th_P_em2->Draw();
    c1.Print(Form("Figs/th_P_em2_%d.pdf", run));
    c1.Print(Form("Figs/th_P_em2_%d.png", run));
    c1.Print(Form("Figs/th_P_em2_%d.root", run));

    TH1D *h_Mmis1 = (TH1D*) file_in.Get("h_Mmis1");
    h_Mmis1->SetLineColor(1);
    h_Mmis1->SetLineWidth(2);
    h_Mmis1->SetTitle("; M_{X}^{2} [GeV^{2}]");
    FormatHist(h_Mmis1);

    TH1D *h_Mmis_corr1 = (TH1D*) file_in.Get("h_Mmis_corr1");
    h_Mmis_corr1->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_corr1->SetLineColor(2);
    h_Mmis_corr1->SetLineWidth(2);
    FormatHist(h_Mmis_corr1);

    TH1D *h_Mmis_corr2 = (TH1D*) file_in.Get("h_Mmis_corr2");
    h_Mmis_corr2->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_corr2->SetLineColor(6);
    h_Mmis_corr2->SetLineWidth(2);
    FormatHist(h_Mmis_corr2);

    TH1D *h_Mmis_corr3 = (TH1D*) file_in.Get("h_Mmis_corr3");
    h_Mmis_corr3->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_corr3->SetLineColor(4);
    h_Mmis_corr3->SetLineWidth(2);
    FormatHist(h_Mmis_corr3);

    TH1D *h_Mmis_AngleFixcorr1 = (TH1D*) file_in.Get("h_Mmis_AngleFixcorr1");
    h_Mmis_AngleFixcorr1->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_AngleFixcorr1->SetLineColor(8);
    h_Mmis_AngleFixcorr1->SetLineWidth(2);
    FormatHist(h_Mmis_AngleFixcorr1);

    TLegend *leg1 = new TLegend(0.5, 0.65, 0.91, 0.9);
    leg1->SetBorderSize(0);
    leg1->AddEntry(h_Mmis1, "Uncorrected");
    leg1->AddEntry(h_Mmis_corr1, "Mom. corrected");
    leg1->AddEntry(h_Mmis_corr2, "Mom & #theta corrected");
    leg1->AddEntry(h_Mmis_corr3, "Mom & #theta & #phi corrected");
    leg1->AddEntry(h_Mmis_AngleFixcorr1, "Mom cor & #theta + #phi Fix");
    h_Mmis_AngleFixcorr1->Draw("hist");
    h_Mmis_corr3->Draw("hist Same");
    h_Mmis_corr2->Draw("hist Same");
    h_Mmis_corr1->Draw("hist Same");
    h_Mmis1->Draw("hist Same");
    leg1->Draw();
    c1.Print(Form("Figs/MmisCorretions_%d.pdf", run));
    c1.Print(Form("Figs/MmisCorretions_%d.png", run));
    c1.Print(Form("Figs/MmisCorretions_%d.root", run));


    TLegend *leg_MM2 = new TLegend(0.65, 0.75, 0.92, 0.92);
    leg_MM2->SetBorderSize(0);
    TH1D *h_Mumu_Corr_AngleFix1 = (TH1D*) file_in.Get("h_Mumu_Corr_AngleFix1");
    h_Mumu_Corr_AngleFix1->SetTitle("; M(#mu^{-}#mu^{+}) [GeV]");
    FormatHist(h_Mumu_Corr_AngleFix1);
    h_Mumu_Corr_AngleFix1->SetLineColor(2);
    h_Mumu_Corr_AngleFix1->Draw("hist");

    TH1D *h_Mumu_Corr_AngleFix2 = (TH1D*) file_in.Get("h_Mumu_Corr_AngleFix2");
    h_Mumu_Corr_AngleFix2->SetTitle("; M(#mu^{-}#mu^{+}) [GeV]");
    FormatHist(h_Mumu_Corr_AngleFix2);
    h_Mumu_Corr_AngleFix2->SetLineColor(4);
    h_Mumu_Corr_AngleFix2->Draw("hist Same");
    leg_MM2->AddEntry(h_Mumu_Corr_AngleFix1, "All events");
    leg_MM2->AddEntry(h_Mumu_Corr_AngleFix2, "After M_{X}^{2} cut");
    leg_MM2->Draw();
    c1.Print(Form("Figs/M_Mumu_AngleFix12_%d.pdf", run));
    c1.Print(Form("Figs/M_Mumu_AngleFix12_%d.png", run));
    c1.Print(Form("Figs/M_Mumu_AngleFix12_%d.root", run));

    TH1D *h_Delta_tM1 = (TH1D*) file_in.Get("h_Delta_tM1");
    h_Delta_tM1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM1->SetLineColor(4);
    h_Delta_tM1->SetLineWidth(2);
    FormatHist(h_Delta_tM1);

    TH1D *h_Delta_tM_Constrian1 = (TH1D*) file_in.Get("h_Delta_tM_Constrian1");
    h_Delta_tM_Constrian1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM_Constrian1->SetLineColor(2);
    h_Delta_tM_Constrian1->SetLineWidth(2);
    FormatHist(h_Delta_tM_Constrian1);

    TH1D *h_Delta_tM_FixMom_Constrain1 = (TH1D*) file_in.Get("h_Delta_tM_FixMom_Constrain1");
    h_Delta_tM_FixMom_Constrain1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM_FixMom_Constrain1->SetLineColor(6);
    h_Delta_tM_FixMom_Constrain1->SetLineWidth(2);
    FormatHist(h_Delta_tM_FixMom_Constrain1);

    TH1D *h_Delta_tM_FixEnergy_Constrain1 = (TH1D*) file_in.Get("h_Delta_tM_FixEnergy_Constrain1");
    h_Delta_tM_FixEnergy_Constrain1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM_FixEnergy_Constrain1->SetLineColor(2);
    h_Delta_tM_FixEnergy_Constrain1->SetLineWidth(2);
    FormatHist(h_Delta_tM_FixEnergy_Constrain1);

    TLegend *leg2 = new TLegend(0.55, 0.7, 0.95, 0.92);
    leg2->SetBorderSize(0);
    leg2->AddEntry(h_Delta_tM1, "Unconstrained");
    leg2->AddEntry(h_Delta_tM_FixMom_Constrain1, "Momentum Constrained");
    leg2->AddEntry(h_Delta_tM_FixEnergy_Constrain1, "Energy Constrained");

    h_Delta_tM1->Draw();
    h_Delta_tM_FixMom_Constrain1->Draw("Same");
    h_Delta_tM_FixEnergy_Constrain1->Draw("Same");
    leg2->Draw();
    c1.Print(Form("Figs/tM_Resolutions_%d.pdf", run));
    c1.Print(Form("Figs/tM_Resolutions_%d.png", run));
    c1.Print(Form("Figs/tM_Resolutions_%d.root", run));


    TH2D *h_Qp2_vs_Q2_Rec1 = (TH2D*)file_in.Get("h_Qp2_vs_Q2_Rec1");
    h_Qp2_vs_Q2_Rec1->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}]");
    FormatHist(h_Qp2_vs_Q2_Rec1);
    h_Qp2_vs_Q2_Rec1->Draw();
    c1.Print(Form("Figs/Qp2_vs_Q2_Rec1_Run_%d.pdf", run));
    c1.Print(Form("Figs/Qp2_vs_Q2_Rec1_Run_%d.png", run));
    c1.Print(Form("Figs/Qp2_vs_Q2_Rec1_Run_%d.root", run));

    TH2D *h_Q2_xB_Rec1 = (TH2D*)file_in.Get("h_Q2_xB_Rec1");
    h_Q2_xB_Rec1->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    FormatHist(h_Q2_xB_Rec1);
    h_Q2_xB_Rec1->Draw();
    c1.Print(Form("Figs/Q2_xB1_Run_%d.pdf", run));
    c1.Print(Form("Figs/Q2_xB1_Run_%d.png", run));
    c1.Print(Form("Figs/Q2_xB1_Run_%d.root", run));
    
    // Plots for \xi vs x binning

    TH2D *h_xB_tM_Rec1 = (TH2D*) file_in.Get("h_xB_tM_Rec1");
    h_xB_tM_Rec1->SetTitle("; -t [GeV^{2}]; x_{B}");
    h_xB_tM_Rec1->SetAxisRange(0., 0.5, "Y");
    h_xB_tM_Rec1->Draw();
    c1.Print(Form("Figs/xB_tM1_Run_%d.pdf", run));
    c1.Print(Form("Figs/xB_tM1_Run_%d.png", run));
    c1.Print(Form("Figs/xB_tM1_Run_%d.root", run));


    const double pol = 0.8; // Assuming 80% beam polarization
    const int n_xi_x_bins = 2;
    const int n_xi_x_Q2bins = 2;

    TH2D * h_xi_xxGPD_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH1D * h_Phi_LH_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D * h_Qp2_Q2_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D * h_tM_xB_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D *h_MC_th_P_em_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; //

    TGraphErrors * gr_Asym_xi_x_bins[n_xi_x_bins][n_xi_x_Q2bins];

    TF1 *f_SinX1 = new TF1("f_SinX1", "[0]*TMath::Sin(x*TMath::DegToRad())", 0., 360.);
    
    for (int i = 0; i < n_xi_x_bins; i++) {
        for (int j = 0; j < n_xi_x_Q2bins; j++) {
            h_xi_xxGPD_xi_x_[i][j] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_xi_x_%d_%d", i, j));
            h_xi_xxGPD_xi_x_[i][j]->SetMarkerColor(i + 2);
            h_Phi_LH_xi_x_[i][j] = (TH1D*) file_in.Get(Form("h_Phi_LH_xi_x_%d_%d", i, j));
            h_Qp2_Q2_xi_x_[i][j] = (TH2D*) file_in.Get(Form("h_Qp2_Q2_xi_x_%d_%d", i, j));
            h_Qp2_Q2_xi_x_[i][j]->SetMarkerColor(i + 1 + 10 * j);
            h_Qp2_Q2_xi_x_[i][j]->SetMarkerStyle(i + 21);
            h_Qp2_Q2_xi_x_[i][j]->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}]");
            h_tM_xB_xi_x_[i][j] = (TH2D*) file_in.Get(Form("h_tM_xB_xi_x_%d_%d", i, j));
            h_tM_xB_xi_x_[i][j]->SetMarkerColor(i + 2);
            h_tM_xB_xi_x_[i][j]->SetTitle("; -t [GeV^{2}]; x_{B}");
            h_MC_th_P_em_xi_x_[i][j] = (TH2D*)file_in.Get(Form("h_MC_th_P_em_xi_x_%d_%d", i, j));
            h_MC_th_P_em_xi_x_[i][j]->SetTitle("; P(e^{-}) [GeV]; #theta (e^{-}) [GeV]");
            h_MC_th_P_em_xi_x_[i][j]->SetMarkerColor(m_cols[i][j]);
            h_MC_th_P_em_xi_x_[i][j]->SetMarkerStyle(21);
            h_MC_th_P_em_xi_x_[i][j]->SetMarkerSize(0.4);
            
            
            double Q2 = h_Qp2_Q2_xi_x_[i][j]->GetMean(1);
            double Qp2 = h_Qp2_Q2_xi_x_[i][j]->GetMean(2);
            double tM = h_tM_xB_xi_x_[i][j]->GetMean(1);
            double xB = h_tM_xB_xi_x_[i][j]->GetMean(2);

            gr_Asym_xi_x_bins[i][j] = new TGraphErrors();
            gr_Asym_xi_x_bins[i][j]->SetMarkerColor(4);
            gr_Asym_xi_x_bins[i][j]->SetMarkerStyle(20);
            gr_Asym_xi_x_bins[i][j]->SetMarkerSize(2);
            gr_Asym_xi_x_bins[i][j]->SetTitle("; #phi_{LH} [deg]; BSA");


            ifstream inp_asym(Form("Dumps/Asym_xi_x_bins_%d_%d_Run_%d.dat", i, j, run));
            ofstream histDump(Form("Dumps/Phi_LH_xi_x_Q2_%d_%d_Run_%d", i, j, run));
            ofstream out_Asum(Form("Dumps/AsymWithErrors_xi_x_bins_%d_%d_Run_%d.dat, run", i, j));

            c1.SetTopMargin(0.1);

            for (int bin = 0; bin < h_Phi_LH_xi_x_[i][j]->GetNbinsX(); bin++) {
                double N_evBin = h_Phi_LH_xi_x_[i][j]->GetBinContent(bin + 1);
                histDump << setw(3) << h_Phi_LH_xi_x_[i][j]->GetBinCenter(bin + 1) << setw(9) << N_evBin << endl;

                double phi, asym;
                inp_asym>> phi;
                inp_asym>> asym;

                double asymErr = (1. / pol) * sqrt((1 - pol * asym * pol * asym) / N_evBin);

                gr_Asym_xi_x_bins[i][j]->SetPoint(bin, phi, asym);
                gr_Asym_xi_x_bins[i][j]->SetPointError(bin, 0, asymErr);
                out_Asum << setw(5) << phi << setw(12) << asym << setw(12) << asymErr << endl;
            }

            gr_Asym_xi_x_bins[i][j]->Draw("AP");
            gr_Asym_xi_x_bins[i][j]->Fit(f_SinX1, "MeV", "", 0., 360.);
            double pol = f_SinX1->GetParameter(0);
            double polErr = f_SinX1->GetParError(0);
            lat1->DrawLatex(0.04, 0.93, Form("Q^{2}=%1.2f GeV^{2}  Q^{'2} = %1.2f GeV^{2}  xB=%1.2f  -t=%1.2f GeV^{2}", Q2, Qp2, xB, tM));
            lat1->DrawLatex(0.25, 0.5 - 0.35*pol/TMath::Abs(pol), Form("pol = %1.1f %% #pm %1.1f %%", 100*pol, 100*polErr));
            c1.Print(Form("Figs/Asym_xi_x_bins_%d_%d_Run_%d.pdf", i, j, run));
            c1.Print(Form("Figs/Asym_xi_x_bins_%d_%d_Run_%d.png", i, j, run));
            c1.Print(Form("Figs/Asym_xi_x_bins_%d_%d_Run_%d.root", i, j, run));

        }
    }


    TLine *line1 = new TLine();
    line1->SetLineColor(2);
    line1->SetLineWidth(2);
    line1->SetLineStyle(9);
    TH2D *h_xi_xxGPD_Rec1 = (TH2D*) file_in.Get("h_xi_xxGPD_Rec1");
    h_xi_xxGPD_Rec1->SetTitle("; x; #xi");
    h_xi_xxGPD_Rec1->SetAxisRange(0., 0.4, "Y");
    h_xi_xxGPD_Rec1->SetAxisRange(-0.2, 0.2, "X");
    h_xi_xxGPD_Rec1->Draw("Scat");
    //h_xi_xxGPD_xi_x_[0][0]->SetMarkerSize(3);
    h_xi_xxGPD_xi_x_[0][0]->SetMarkerStyle(3);
    h_xi_xxGPD_xi_x_[0][0]->Draw("Same scat");
    h_xi_xxGPD_xi_x_[1][0]->SetMarkerStyle(20);
    //h_xi_xxGPD_xi_x_[1][0]->SetMarkerSize(1);
    h_xi_xxGPD_xi_x_[1][0]->Draw("Same scat");
    line1->DrawLine(0., 0., 0.2, 0.2);
    line1->DrawLine(0., 0., -0.2, 0.2);
    c1.Print(Form("Figs/Xi_x_bins_on_xi_x_Run_%d.pdf", run));
    c1.Print(Form("Figs/Xi_x_bins_on_xi_x_Run_%d.png", run));
    c1.Print(Form("Figs/Xi_x_bins_on_xi_x_Run_%d.root", run));


    h_Qp2_Q2_xi_x_[0][0]->Draw("Scat");
    h_Qp2_Q2_xi_x_[0][1]->SetMarkerColor(9);
    h_Qp2_Q2_xi_x_[0][1]->Draw("Same Scat");
    h_Qp2_Q2_xi_x_[1][0]->Draw("Same Scat");
    h_Qp2_Q2_xi_x_[1][1]->Draw("Same Scat");
    c1.Print(Form("Figs/Qp2_Q2_With_xi_x_bins_Run_%d.pdf", run));
    c1.Print(Form("Figs/Qp2_Q2_With_xi_x_bins_Run_%d.png", run));
    c1.Print(Form("Figs/Qp2_Q2_With_xi_x_bins_Run_%d.root", run));

    h_tM_xB_xi_x_[0][0]->Draw("Scat");
    h_tM_xB_xi_x_[1][0]->Draw("Same Scat");
    c1.Print(Form("Figs/xB_tM_xi_x_bins_Run_%d.pdf", run));
    c1.Print(Form("Figs/xB_tM_xi_x_bins_Run_%d.png", run));
    c1.Print(Form("Figs/xB_tM_xi_x_bins_Run_%d.root", run));

    
    TH2D *h_MC_th_P_em_FS1 = (TH2D*)file_in.Get("h_MC_th_P_em_FS1");
    h_MC_th_P_em_FS1->SetTitle("; P(e^{-}) [GeV]; #theta (e^{-}) [GeV]");
    h_MC_th_P_em_FS1->SetAxisRange(0., 35., "Y");
    h_MC_th_P_em_FS1->SetAxisRange(0., 8., "X");
    h_MC_th_P_em_FS1->Draw("Scat");
    h_MC_th_P_em_xi_x_[0][0]->Draw("Same Scat");
    h_MC_th_P_em_xi_x_[0][1]->Draw("Same Scat");
    h_MC_th_P_em_xi_x_[1][0]->Draw("Same Scat");
    h_MC_th_P_em_xi_x_[1][1]->Draw("Same Scat");
    lat1->SetTextColor(m_cols[0][0]);
    lat1->DrawLatex(0.75, 0.86, "Bin(0,0)");
    lat1->SetTextColor(m_cols[0][1]);
    lat1->DrawLatex(0.75, 0.81, "Bin(0,1)");
    lat1->SetTextColor(m_cols[1][0]);
    lat1->DrawLatex(0.75, 0.76, "Bin(1,0)");
    lat1->SetTextColor(m_cols[1][1]);
    lat1->DrawLatex(0.75, 0.71, "Bin(1,1)");
    
    c1.Print(Form("Figs/th_P_em_xi_x_bins_%d.pdf", run));
    c1.Print(Form("Figs/th_P_em_xi_x_bins_%d.png", run));
    c1.Print(Form("Figs/th_P_em_xi_x_bins_%d.root", run));
    
    return 0;
}

template <typename T>
void FormatHist(T h) {

    //TH2D *h = new TH2D("dsd", "", 100, 0, 100, 200, 0, 100);

    if (h->GetEntries() < 25) {
        h->SetMinimum(-0.1);
    }

    h->SetTitleSize(0.05, "X");
    h->SetTitleSize(0.05, "Y");
    h->SetLabelSize(0.05, "X");
    h->SetLabelSize(0.05, "Y");
}