/* 
 * File:   DrawPlots.cc
 * Author: rafopar
 *
 * Created on June 4, 2022, 10:09 PM
 */

#include <TH1D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <cstdlib>

using namespace std;

/*
 * 
 */
int DrawPlots(int run) {

    const int nQp2Bins_Bin1 = 4;
    const int nQ2Bins_Bin1 = 4;
    const int nQp2Bins_Bin2 = 4;
    const int nQp2Bins_Bin3 = 3;
    const double pol = 0.8; // Assuming 80% beam polarization

    /**
     * Defining Bin 1
     */

    const double t_max_Bin1 = 0.4;
    const double t_min_Bin1 = 0.1;
    const double xB_max_Bin1 = 0.22;
    const double xB_min_Bin1 = 0.12;
    const double Qp2_edges[nQp2Bins_Bin1 + 1] = {0.8, 1.6, 2.4, 3.2, 4.};
    const double Q2_edges[nQ2Bins_Bin1 + 1] = {1., 1.6, 2.4, 3.2, 4.};
    const double Qp2_max_bin1 = 3.;
    const double Qp2_min_bin1 = 2.;
    const double Q2_max_bin1 = 3.;
    const double Q2_min_bin1 = 2.;

    /**
     * Defining Bin 2
     */
    const double t_max_Bin2 = 0.4;
    const double t_min_Bin2 = 0.05;
    const double xB_max_Bin2 = 0.01;
    const double xB_min_Bin2 = 0.05;

    const double Q2_max_bin2 = 0.84;
    const double Q2_min_bin2 = 0.15;
    const double Qp2_edges_bin2[nQp2Bins_Bin2 + 1] = {1., 2., 3., 4., 6.};

    /**
     * Defining Bin 3
     */
    const double t_max_Bin3 = 0.4;
    const double t_min_Bin3 = 0.05;
    const double xB_max_Bin3 = 0.05;
    const double xB_min_Bin3 = 0.12;

    const double Q2_max_bin3 = 2.15;
    const double Q2_min_bin3 = 0.65;
    const double Qp2_edges_bin3[nQp2Bins_Bin3 + 1] = {1., 2., 3., 5.};

    TLine *line1 = new TLine();

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);

    TFile file_in(Form("AnaGrape_%d.root", run));

    c1->SetGridy();
    c1->SetGridx();
    c1->SetLogz();
    TH2D *h_Qp2_vs_Q2_2 = (TH2D*) file_in.Get("h_Qp2_vs_Q2_2");
    h_Qp2_vs_Q2_2->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2}[GeV^{2}]");
    h_Qp2_vs_Q2_2->Draw("colz");
    c1->Print(Form("Figs/Qp2_vs_Q2_2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_2_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_2_Run_%d.root", run));

    TH2D *h_Q2_xB2 = (TH2D*) file_in.Get("h_Q2_xB2");
    h_Q2_xB2->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2->Draw();
    c1->Print(Form("Figs/Q2_xB2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Q2_xB2_Run_%d.png", run));
    c1->Print(Form("Figs/Q2_xB2_Run_%d.root", run));

    c1->SetLogz(0);
    TH2D *h_th_P_p2 = (TH2D*) file_in.Get("h_th_P_p2");
    h_th_P_p2->SetTitle("; P [GeV]; #theta [deg]");
    h_th_P_p2->Draw("colz");
    c1->Print(Form("Figs/th_P_p2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_p2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_p2_Run_%d.root", run));

    TH2D *h_th_P_p3 = (TH2D*) file_in.Get("h_th_P_p3");
    h_th_P_p3->SetTitle("; P [GeV]; #theta [deg]");
    h_th_P_p3->Draw("colz");
    c1->Print(Form("Figs/th_P_p3_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_p3_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_p3_Run_%d.root", run));

    TH2D *h_th_P_em2 = (TH2D*) file_in.Get("h_th_P_em2");
    h_th_P_em2->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-}) [deg]");
    h_th_P_em2->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_em2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_em2_Run_%d.root", run));

    TH2D *h_th_P_mum2 = (TH2D*) file_in.Get("h_th_P_mum2");
    h_th_P_mum2->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) [deg]");
    h_th_P_mum2->Draw("colz");
    c1->Print(Form("Figs/th_P_mum2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_mum2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_mum2_Run_%d.root", run));

    TH2D *h_th_P_mup2 = (TH2D*) file_in.Get("h_th_P_mup2");
    h_th_P_mup2->SetTitle("; P(#mu^{+}) [GeV]; #theta(#mu^{+}) [deg]");
    h_th_P_mup2->Draw("colz");
    c1->Print(Form("Figs/th_P_mup2_Run_%d.pdf", run));
    c1->Print(Form("Figs/th_P_mup2_Run_%d.png", run));
    c1->Print(Form("Figs/th_P_mup2_Run_%d.root", run));

    line1->SetLineColor(95);
    line1->SetLineWidth(5);
    TH2D *h_xB_tM2 = (TH2D*) file_in.Get("h_xB_tM2");
    h_xB_tM2->SetTitle("; -t [GeV^{2}]; x_{B}");
    h_xB_tM2->SetAxisRange(0., 0.7, "Y");
    h_xB_tM2->Draw("colz");
    line1->DrawLine(t_min_Bin1, xB_min_Bin1, t_max_Bin1, xB_min_Bin1);
    line1->DrawLine(t_max_Bin1, xB_min_Bin1, t_max_Bin1, xB_max_Bin1);
    line1->DrawLine(t_max_Bin1, xB_max_Bin1, t_min_Bin1, xB_max_Bin1);
    line1->DrawLine(t_min_Bin1, xB_max_Bin1, t_min_Bin1, xB_min_Bin1);
    line1->DrawLine(t_min_Bin2, xB_min_Bin2, t_max_Bin2, xB_min_Bin2);
    line1->DrawLine(t_max_Bin2, xB_min_Bin2, t_max_Bin2, xB_max_Bin2);
    line1->DrawLine(t_max_Bin2, xB_max_Bin2, t_min_Bin2, xB_max_Bin2);
    line1->DrawLine(t_min_Bin2, xB_max_Bin2, t_min_Bin2, xB_min_Bin2);
    line1->DrawLine(t_min_Bin3, xB_min_Bin3, t_max_Bin3, xB_min_Bin3);
    line1->DrawLine(t_max_Bin3, xB_min_Bin3, t_max_Bin3, xB_max_Bin3);
    line1->DrawLine(t_max_Bin3, xB_max_Bin3, t_min_Bin3, xB_max_Bin3);
    line1->DrawLine(t_min_Bin3, xB_max_Bin3, t_min_Bin3, xB_min_Bin3);

    c1->Print(Form("Figs/xB_tM2_Run_%d.pdf", run));
    c1->Print(Form("Figs/xB_tM2_Run_%d.png", run));
    c1->Print(Form("Figs/xB_tM2_Run_%d.root", run));

    TH2D *h_Qp2_vs_Q2_Bin1 = (TH2D*) file_in.Get("h_Qp2_vs_Q2_3");
    h_Qp2_vs_Q2_Bin1->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2}[GeV^{2}]");
    h_Qp2_vs_Q2_Bin1->Draw("colz");
    line1->DrawLine(Q2_min_bin1, Qp2_edges[0], Q2_min_bin1, Qp2_edges[nQp2Bins_Bin1]);
    line1->DrawLine(Q2_max_bin1, Qp2_edges[0], Q2_max_bin1, Qp2_edges[nQp2Bins_Bin1]);
    for (int i = 0; i < nQp2Bins_Bin1 + 1; i++) {
        line1->DrawLine(Q2_min_bin1, Qp2_edges[i], Q2_max_bin1, Qp2_edges[i]);
    }
    line1->SetLineColor(6);
    line1->DrawLine(Q2_edges[0], Qp2_min_bin1, Q2_edges[nQ2Bins_Bin1], Qp2_min_bin1);
    line1->DrawLine(Q2_edges[0], Qp2_max_bin1, Q2_edges[nQ2Bins_Bin1], Qp2_max_bin1);
    for (int i = 0; i < nQ2Bins_Bin1 + 1; i++) {
        line1->DrawLine(Q2_edges[i], Qp2_min_bin1, Q2_edges[i], Qp2_max_bin1);
    }

    c1->Print(Form("Figs/Qp2_vs_Q2_Bin1_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin1_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin1_Run_%d.root", run));
    // const double Qp2_edges[nQp2Bins_Bin1 + 1] = {0.8, 1.6, 2.4, 3.2, 4.};
    // Q2_min_bin1

    TH2D *h_Qp2_vs_Q2_Bin2 = (TH2D*) file_in.Get("h_Qp2_vs_Q2_4");
    h_Qp2_vs_Q2_Bin2->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2}[GeV^{2}]");
    h_Qp2_vs_Q2_Bin2->Draw("colz");
    line1->SetLineColor(95);
    line1->DrawLine(Q2_min_bin2, Qp2_edges_bin2[0], Q2_min_bin2, Qp2_edges_bin2[nQp2Bins_Bin2]);
    line1->DrawLine(Q2_max_bin2, Qp2_edges_bin2[0], Q2_max_bin2, Qp2_edges_bin2[nQp2Bins_Bin2]);
    for (int i = 0; i < nQp2Bins_Bin2 + 1; i++) {
        line1->DrawLine(Q2_min_bin2, Qp2_edges_bin2[i], Q2_max_bin2, Qp2_edges_bin2[i]);
    }
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin2_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin2_Run_%d.root", run));

    TH2D *h_Qp2_vs_Q2_Bin3 = (TH2D*) file_in.Get("h_Qp2_vs_Q2_5");
    h_Qp2_vs_Q2_Bin3->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2}[GeV^{2}]");
    h_Qp2_vs_Q2_Bin3->Draw("colz");
    line1->DrawLine(Q2_min_bin3, Qp2_edges_bin3[0], Q2_min_bin3, Qp2_edges_bin3[nQp2Bins_Bin3]);
    line1->DrawLine(Q2_max_bin3, Qp2_edges_bin3[0], Q2_max_bin3, Qp2_edges_bin3[nQp2Bins_Bin3]);
    for (int i = 0; i < nQp2Bins_Bin3 + 1; i++) {
        line1->DrawLine(Q2_min_bin3, Qp2_edges_bin3[i], Q2_max_bin3, Qp2_edges_bin3[i]);
    }
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin3_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin3_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_Bin3_Run_%d.root", run));


    TH2D *h_xi_xxGPD2 = (TH2D*) file_in.Get("h_xi_xxGPD2");
    h_xi_xxGPD2->SetTitle("; x; #xi");
    h_xi_xxGPD2->SetAxisRange(0., 0.4, "Y");
    h_xi_xxGPD2->SetAxisRange(-0.4, 0.4, "X");
    h_xi_xxGPD2->Draw("Scat");

    TH2D * h_xi_xxGPD2_[nQp2Bins_Bin1];

    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.14);

    for (int i = 0; i < nQp2Bins_Bin1; i++) {
        h_xi_xxGPD2_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD2_%d", i));
        h_xi_xxGPD2_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD2_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD2_[i]->Draw("Same Scat");
    }
    c1->Print(Form("Figs/xi_vs_x_FixQ2_RUn_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_FixQ2_RUn_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_FixQ2_RUn_%d.root", run));

    TH2D * h_xi_xxGPD3_[nQp2Bins_Bin1];

    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < nQp2Bins_Bin1; i++) {
        h_xi_xxGPD3_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD3_%d", i));
        h_xi_xxGPD3_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD3_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD3_[i]->Draw("Same Scat");
    }
    c1->Print(Form("Figs/xi_vs_x_FixQp2_RUn_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_FixQp2_RUn_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_FixQp2_RUn_%d.root", run));


    // Kinematic Bin2
    TH2D * h_xi_xxGPD_bin2_[nQp2Bins_Bin2];
    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < nQp2Bins_Bin2; i++) {
        h_xi_xxGPD_bin2_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_bin2_%d", i));
        h_xi_xxGPD_bin2_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD_bin2_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD_bin2_[i]->SetMarkerSize(2);
        h_xi_xxGPD_bin2_[i]->Draw("Same Scat");
    }
    c1->Print(Form("Figs/xi_vs_x_KinBin2_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin2_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin2_Run_%d.root", run));


    // Kinematic Bin3
    TH2D * h_xi_xxGPD_bin3_[nQp2Bins_Bin3];
    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < nQp2Bins_Bin3; i++) {
        h_xi_xxGPD_bin3_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_bin3_%d", i));
        h_xi_xxGPD_bin3_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD_bin3_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD_bin3_[i]->Draw("Same Scat");
    }
    c1->Print(Form("Figs/xi_vs_x_KinBin3_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin3_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin3_Run_%d.root", run));



    // ******* New kinematic bins, that we chose with Stepan in July 2024 for the 22 GeV Discusson meeting.

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);
    lat1->SetTextColor(4);
    lat1->SetTextSize(0.04);

    //TLine *line1 = new TLine();
    line1->SetLineColor(95);
    line1->SetLineWidth(3);

    const int n_Q2bins = 4;
    const int n_Qp2bins = 4;

    const int Qp2_Max = 3.;
    const int Qp2_Min = 2.;
    const int Q2_Max = 2;
    const int Q2_Min = 1;
    const double Q2_edges_new[n_Q2bins + 1] = {1., 1.3, 1.7, 2.5, 4.};
    const double Qp2_edges_new[n_Q2bins + 1] = {2., 2.3, 2.9, 3.7, 5.};
    const double tMcut = 0.4;

    h_xB_tM2->Draw("colz");
    line1->DrawLine(tMcut, 0., tMcut, 0.6);
    c1->Print(Form("Figs/xB_tM2new_Run_%d.pdf", run));
    c1->Print(Form("Figs/xB_tM2new_Run_%d.png", run));
    c1->Print(Form("Figs/xB_tM2new_Run_%d.root", run));

    TH2D *h_Qp2_vs_Q2_6 = (TH2D*) file_in.Get("h_Qp2_vs_Q2_6");
    h_Qp2_vs_Q2_6->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2}[GeV^{2}]");
    h_Qp2_vs_Q2_6->Draw("colz");

    line1->DrawLine(Q2_edges_new[0], Qp2_Max, Q2_edges_new[n_Q2bins], Qp2_Max);
    line1->DrawLine(Q2_edges_new[0], Qp2_Min, Q2_edges_new[n_Q2bins], Qp2_Min);
    for (int i = 0; i < n_Q2bins + 1; i++) {
        line1->DrawLine(Q2_edges_new[i], Qp2_Min, Q2_edges_new[i], Qp2_Max);
    }

    c1->Print(Form("Figs/Qp2_vs_Q2_6_Q2Scan_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Q2Scan_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Q2Scan_Run_%d.root", run));

    h_Qp2_vs_Q2_6->Draw("colz");
    line1->SetLineColor(6);
    line1->SetLineWidth(3);
    line1->DrawLine(Q2_Min, Qp2_edges_new[0], Q2_Min, Qp2_edges_new[n_Qp2bins]);
    line1->DrawLine(Q2_Max, Qp2_edges_new[0], Q2_Max, Qp2_edges_new[n_Qp2bins]);
    for (int i = 0; i < n_Q2bins + 1; i++) {
        line1->DrawLine(Q2_Min, Qp2_edges_new[i], Q2_Max, Qp2_edges_new[i]);
    }
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Qp2Scan_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Qp2Scan_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Qp2Scan_Run_%d.root", run));


    TH2D * h_xi_xxGPD_newQ2Scan_[n_Q2bins];
    TH1D * h_Phi_LH_newQ2Scan_[n_Q2bins];
    TH2D * h_Qp2_Q2_newQ2Scan_[n_Q2bins];
    TH2D * h_tM_xB_newQ2Scan_[n_Q2bins];
    TGraphErrors * gr_A_LH_newQ2Scan_[n_Q2bins];

    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < n_Q2bins; i++) {
        h_xi_xxGPD_newQ2Scan_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_newQ2Scan_%d", i));
        h_xi_xxGPD_newQ2Scan_[i]->SetMarkerColor(i + 2);
        lat1->SetTextColor(i + 2);
        h_xi_xxGPD_newQ2Scan_[i]->Draw("Same scat");
        lat1->DrawLatex(0.65, 0.4 - 0.05 * i, Form("Q^{2}#in (%1.1f-%1.1f)GeV", Q2_edges_new[i], Q2_edges_new[i + 1]));
    }
    c1->Print(Form("Figs/xi_vs_x_New_Q2Scan_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_New_Q2Scan_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_New_Q2Scan_Run_%d.root", run));

    ofstream out_Q2Scan(Form("Dumps/Q2_Scan_Kinematics_%d.dat", run));
    out_Q2Scan << "bin" << setw(10) << "Q2" << setw(12) << "Qp2" << setw(12) << "-t" << setw(12) << "xB" << endl;

    lat1->SetTextColor(4);
    for (int i = 0; i < n_Q2bins; i++) {
        h_Phi_LH_newQ2Scan_[i] = (TH1D*) file_in.Get(Form("h_Phi_LH_newQ2Scan_%d", i));
        h_Phi_LH_newQ2Scan_[i]->SetMinimum(0);
        h_Phi_LH_newQ2Scan_[i]->SetTitle("; #phi_{LH} [deg]");
        h_Phi_LH_newQ2Scan_[i]->Draw();
        c1->Print(Form("Figs/Phi_LH_Q2Scan_Bin_%d_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Phi_LH_Q2Scan_Bin_%d_Run_%d.png", i, run));
        c1->Print(Form("Figs/Phi_LH_Q2Scan_Bin_%d_Run_%d.root", i, run));
        ofstream histDump(Form("Dumps/Phi_LH_Q2Scan_Bin_%d_Run_%d.dat", i, run));
        ifstream inp_asym(Form("Dumps/Asym_Q2Scan_Bin%d.dat", i));

        gr_A_LH_newQ2Scan_[i] = new TGraphErrors();
        gr_A_LH_newQ2Scan_[i]->SetMarkerColor(4);
        gr_A_LH_newQ2Scan_[i]->SetMarkerStyle(20);
        gr_A_LH_newQ2Scan_[i]->SetMarkerSize(2);
        gr_A_LH_newQ2Scan_[i]->SetTitle("; #phi_{LH} [deg]; BSA");

        histDump << setw(3) << "Phi" << setw(13) << "# of events" << endl;
        for (int bin = 0; bin < h_Phi_LH_newQ2Scan_[i]->GetNbinsX(); bin++) {
            double N_evBin = h_Phi_LH_newQ2Scan_[i]->GetBinContent(bin + 1);
            histDump << setw(3) << h_Phi_LH_newQ2Scan_[i]->GetBinCenter(bin + 1) << setw(9) << N_evBin << endl;

            double phi, asym;
            inp_asym>> phi;
            inp_asym>> asym;

            double asymErr = (1. / pol) * sqrt((1 - pol * asym * pol * asym) / N_evBin);

            gr_A_LH_newQ2Scan_[i]->SetPoint(bin, phi, asym);
            gr_A_LH_newQ2Scan_[i]->SetPointError(bin, 0, asymErr);
        }

        c1->Clear();

        h_Qp2_Q2_newQ2Scan_[i] = (TH2D*) file_in.Get(Form("h_Qp2_Q2_newQ2Scan_%d", i));
        double Q2 = h_Qp2_Q2_newQ2Scan_[i]->GetMean(1);
        double Qp2 = h_Qp2_Q2_newQ2Scan_[i]->GetMean(2);
        h_tM_xB_newQ2Scan_[i] = (TH2D*) file_in.Get(Form("h_tM_xB_newQ2Scan_%d", i));
        double tM = h_tM_xB_newQ2Scan_[i]->GetMean(1);
        double xB = h_tM_xB_newQ2Scan_[i]->GetMean(2);
        out_Q2Scan << i << setw(12) << Q2 << setw(12) << Qp2 << setw(12) << tM << setw(12) << xB << endl;

        gr_A_LH_newQ2Scan_[i]->Draw("AP");
        lat1->DrawLatex(0.04, 0.93, Form("Q^{2}=%1.2f GeV^{2}  Q^{'2} = %1.2f GeV^{2}  xB=%1.2f  -t=%1.2f GeV^{2}", Q2, Qp2, xB, tM));

        c1->Print(Form("Figs/Asym_Q2Scan_Bin_%d_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Asym_Q2Scan_Bin_%d_Run_%d.png", i, run));
        c1->Print(Form("Figs/Asym_Q2Scan_Bin_%d_Run_%d.root", i, run));
    }

    TH2D * h_xi_xxGPD_newQp2Scan_[n_Q2bins];
    TH1D * h_Phi_LH_newQp2Scan_[n_Q2bins];
    TH2D * h_Qp2_Q2_newQp2Scan_[n_Q2bins];
    TH2D * h_tM_xB_newQp2Scan_[n_Q2bins];
    TGraphErrors * gr_A_LH_newQp2Scan_[n_Q2bins];

    ofstream out_Qp2Scan(Form("Dumps/Qp2_Scan_Kinematics_%d.dat", run));




    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < n_Qp2bins; i++) {
        h_xi_xxGPD_newQp2Scan_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_newQp2Scan_%d", i));
        h_xi_xxGPD_newQp2Scan_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD_newQp2Scan_[i]->Draw("Same scat");
        lat1->SetTextColor(i + 2);
        lat1->DrawLatex(0.65, 0.4 - 0.05 * i, Form("Q^{'2}#in (%1.1f-%1.1f)GeV", Qp2_edges_new[i], Qp2_edges_new[i + 1]));
    }
    c1->Print(Form("Figs/xi_vs_x_New_Qp2Scan_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_New_Qp2Scan_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_New_Qp2Scan_Run_%d.root", run));

    lat1->SetTextColor(4);
    out_Qp2Scan << "bin" << setw(10) << "Q2" << setw(12) << "Qp2" << setw(12) << "-t" << setw(12) << "xB" << endl;
    for (int i = 0; i < n_Qp2bins; i++) {
        h_Phi_LH_newQp2Scan_[i] = (TH1D*) file_in.Get(Form("h_Phi_LH_newQp2Scan_%d", i));
        h_Phi_LH_newQp2Scan_[i]->SetMinimum(0);
        h_Phi_LH_newQp2Scan_[i]->SetTitle("; #phi_{LH} [deg]");
        h_Phi_LH_newQp2Scan_[i]->Draw();
        c1->Print(Form("Figs/Phi_LH_Qp2Scan_Bin_%d_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Phi_LH_Qp2Scan_Bin_%d_Run_%d.png", i, run));
        c1->Print(Form("Figs/Phi_LH_Qp2Scan_Bin_%d_Run_%d.root", i, run));

        ofstream histDump(Form("Dumps/Phi_LH_Qp2Scan_Bin_%d_Run_%d.dat", i, run));
        ifstream inp_asym(Form("Dumps/Asym_Qp2Scan_Bin%d.dat", i));

        gr_A_LH_newQp2Scan_[i] = new TGraphErrors();
        gr_A_LH_newQp2Scan_[i]->SetMarkerColor(4);
        gr_A_LH_newQp2Scan_[i]->SetMarkerStyle(20);
        gr_A_LH_newQp2Scan_[i]->SetMarkerSize(2);
        gr_A_LH_newQp2Scan_[i]->SetTitle("; #phi_{LH} [deg]; BSA");

        histDump << setw(3) << "Phi" << setw(13) << "# of events" << endl;
        for (int bin = 0; bin < h_Phi_LH_newQp2Scan_[i]->GetNbinsX(); bin++) {
            double N_evBin = h_Phi_LH_newQp2Scan_[i]->GetBinContent(bin + 1);
            histDump << setw(3) << h_Phi_LH_newQp2Scan_[i]->GetBinCenter(bin + 1) << setw(9) << N_evBin << endl;

            double phi, asym;
            inp_asym>> phi;
            inp_asym>> asym;

            double asymErr = (1. / pol) * sqrt((1 - pol * asym * pol * asym) / N_evBin);

            gr_A_LH_newQ2Scan_[i]->SetPoint(bin, phi, asym);
            gr_A_LH_newQ2Scan_[i]->SetPointError(bin, 0, asymErr);
        }

        c1->Clear();

        h_Qp2_Q2_newQp2Scan_[i] = (TH2D*) file_in.Get(Form("h_Qp2_Q2_newQp2Scan_%d", i));
        double Q2 = h_Qp2_Q2_newQp2Scan_[i]->GetMean(1);
        double Qp2 = h_Qp2_Q2_newQp2Scan_[i]->GetMean(2);
        h_tM_xB_newQp2Scan_[i] = (TH2D*) file_in.Get(Form("h_tM_xB_newQp2Scan_%d", i));
        double tM = h_tM_xB_newQp2Scan_[i]->GetMean(1);
        double xB = h_tM_xB_newQp2Scan_[i]->GetMean(2);
        out_Qp2Scan << i << setw(12) << Q2 << setw(12) << Qp2 << setw(12) << tM << setw(12) << xB << endl;

        gr_A_LH_newQ2Scan_[i]->Draw("AP");
        lat1->DrawLatex(0.04, 0.93, Form("Q^{2}=%1.2f GeV^{2}  Q^{'2} = %1.2f GeV^{2}  xB=%1.2f  -t=%1.2f GeV^{2}", Q2, Qp2, xB, tM));

        c1->Print(Form("Figs/Asym_Qp2Scan_Bin_%d_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Asym_Qp2Scan_Bin_%d_Run_%d.png", i, run));
        c1->Print(Form("Figs/Asym_Qp2Scan_Bin_%d_Run_%d.root", i, run));

    }


    return 0;
}