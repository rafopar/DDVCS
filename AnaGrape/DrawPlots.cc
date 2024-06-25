/* 
 * File:   DrawPlots.cc
 * Author: rafopar
 *
 * Created on June 4, 2022, 10:09 PM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
int DrawPlots() {

    const int run = 2;

    const int nQp2Bins_Bin1 = 4;
    const int nQ2Bins_Bin1 = 4;
    const int nQp2Bins_Bin2 = 4;
    const int nQp2Bins_Bin3 = 3;


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
    h_xi_xxGPD2->Draw();

    TH2D * h_xi_xxGPD2_[nQp2Bins_Bin1];

    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.14);

    for (int i = 0; i < nQp2Bins_Bin1; i++) {
        h_xi_xxGPD2_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD2_%d", i));
        h_xi_xxGPD2_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD2_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD2_[i]->Draw("Same");
    }
    c1->Print(Form("Figs/xi_vs_x_FixQ2_RUn_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_FixQ2_RUn_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_FixQ2_RUn_%d.root", run));

    TH2D * h_xi_xxGPD3_[nQp2Bins_Bin1];

    h_xi_xxGPD2->Draw();
    for (int i = 0; i < nQp2Bins_Bin1; i++) {
        h_xi_xxGPD3_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD3_%d", i));
        h_xi_xxGPD3_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD3_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD3_[i]->Draw("Same");
    }
    c1->Print(Form("Figs/xi_vs_x_FixQp2_RUn_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_FixQp2_RUn_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_FixQp2_RUn_%d.root", run));


    // Kinematic Bin2
    TH2D * h_xi_xxGPD_bin2_[nQp2Bins_Bin2];
    h_xi_xxGPD2->Draw();
    for (int i = 0; i < nQp2Bins_Bin2; i++) {
        h_xi_xxGPD_bin2_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_bin2_%d", i));
        h_xi_xxGPD_bin2_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD_bin2_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD_bin2_[i]->Draw("Same");
    }
    c1->Print(Form("Figs/xi_vs_x_KinBin2_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin2_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin2_Run_%d.root", run));


    // Kinematic Bin3
    TH2D * h_xi_xxGPD_bin3_[nQp2Bins_Bin3];
    h_xi_xxGPD2->Draw();
    for (int i = 0; i < nQp2Bins_Bin3; i++) {
        h_xi_xxGPD_bin3_[i] = (TH2D*) file_in.Get(Form("h_xi_xxGPD_bin3_%d", i));
        h_xi_xxGPD_bin3_[i]->SetTitle("; x; #xi");
        h_xi_xxGPD_bin3_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD_bin3_[i]->Draw("Same");
    }
    c1->Print(Form("Figs/xi_vs_x_KinBin3_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin3_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_KinBin3_Run_%d.root", run));

    return 0;
}