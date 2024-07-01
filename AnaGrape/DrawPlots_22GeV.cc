/* 
 * File:   DrawPlots_22GeV.cc
 * Author: rafopar
 *
 * Created on July 8, 2022, 10:26 AM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
int DrawPlots_22GeV( int run ) {

    gStyle->SetOptStat(0);

    //const int run = 16;
    int run_11GeV = 2;

    TF1 *f_Q2_FixW = new TF1("f_Q2_FixW", "([0]*[0] - 0.9383*0.9383)*x/(1.-x)", 0., 1.);
    f_Q2_FixW->SetParameter(0, 2.);

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.02);

    TFile *file_in = new TFile(Form("AnaGrape_22GeV_%d.root", run));

    TH2D *h_xB_tM2 = (TH2D*) file_in->Get("h_xB_tM2");
    h_xB_tM2->SetTitle("; -t [GeV^{2}]; x_{B}");
    h_xB_tM2->SetAxisRange(0., 0.7, "Y");
    h_xB_tM2->Draw("colz");
    c1->Print(Form("Figs/xB_tM_22GeV_Run%d.pdf", run));
    c1->Print(Form("Figs/xB_tM_22GeV_Run%d.png", run));
    c1->Print(Form("Figs/xB_tM_22GeV_Run%d.root", run));

    TH2D *h_Qp2_vs_Q2_2 = (TH2D*) file_in->Get("h_Qp2_vs_Q2_2");
    h_Qp2_vs_Q2_2->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}] ");
    h_Qp2_vs_Q2_2->SetAxisRange(2.5, 8.5, "X");
    h_Qp2_vs_Q2_2->Draw("colz");
    c1->Print(Form("Figs/Qp2_Q2_2_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Q2_2_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Q2_2_22GeV_Run_%d.root", run));

    TH2D *h_Qp2_vs_Q2_3 = (TH2D*) file_in->Get("h_Qp2_vs_Q2_3");
    h_Qp2_vs_Q2_3->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}] ");
    //h_Qp2_vs_Q2_3->SetAxisRange(2.5, 8.5, "X");
    h_Qp2_vs_Q2_3->Draw("colz");
    c1->Print(Form("Figs/Qp2_Q2_3_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Q2_3_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Q2_3_22GeV_Run_%d.root", run));

    TH2D *h_Qp2_vs_Q2_6 = (TH2D*) file_in->Get("h_Qp2_vs_Q2_6");
    h_Qp2_vs_Q2_6->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}] ");
    //h_Qp2_vs_Q2_6->SetAxisRange(2.5, 8.5, "X");
    h_Qp2_vs_Q2_6->Draw("colz");
    c1->Print(Form("Figs/Qp2_Q2_6_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Q2_6_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Q2_6_22GeV_Run_%d.root", run));



    const int nQp2Bins = 8;
    const int nQ2Bins = 8;

    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.12);
    TH2D *h_xi_xxGPD2 = (TH2D*) file_in->Get("h_xi_xxGPD2");
    h_xi_xxGPD2->SetAxisRange(0, 0.4, "Y");
    h_xi_xxGPD2->SetAxisRange(-0.4, 0.4, "X");
    h_xi_xxGPD2->SetTitle("; x; #xi");
    h_xi_xxGPD2->Draw();

    TH2D * h_xi_xxGPD2_[nQp2Bins];
    TH2D * h_xi_xxGPD3_[nQ2Bins];
    TH2D * h_xi_xxGPD4_[nQp2Bins];
    TH2D * h_xi_xxGPD5_[nQp2Bins];

    TLegend *leg1 = new TLegend(0.65, 0.6, 0.92, 0.92);
    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextSize(0.025);

    for (int i = 0; i < nQp2Bins; i++) {
        h_xi_xxGPD2_[i] = (TH2D*) file_in->Get(Form("h_xi_xxGPD2_%d", i));
        h_xi_xxGPD2_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD2_[i]->Draw("Same");
        //leg1->AddEntry(h_xi_xxGPD2_[i], Form("%1.2f < Q^{2}[GeV^{2}] < %1.2f", i+1., i+2.) );
        lat1->SetTextColor(i + 2);
        lat1->DrawLatex(0.72, 0.12 + i * 0.03, Form("%1.2f < Q^{'2}[GeV^{2}] < %1.2f", i + 1., i + 2.));
    }
    lat1->SetTextColor(1);
    lat1->DrawLatex(0.17, 0.2, "3.2 < Q^{2} [GeV] < 4.2");
    c1->Print(Form("Figs/Qp2_Scan1_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Scan1_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Scan1_Run_%d.root", run));

    h_xi_xxGPD2->Draw();
    for (int i = 0; i < nQp2Bins; i++) {
        h_xi_xxGPD4_[i] = (TH2D*) file_in->Get(Form("h_xi_xxGPD4_%d", i));
        h_xi_xxGPD4_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD4_[i]->Draw("Same");
        //leg1->AddEntry(h_xi_xxGPD2_[i], Form("%1.2f < Q^{2}[GeV^{2}] < %1.2f", i+1., i+2.) );
        lat1->SetTextColor(i + 2);
        lat1->DrawLatex(0.72, 0.12 + i * 0.03, Form("%1.2f < Q^{'2}[GeV^{2}] < %1.2f", i + 1., i + 2.));
    }
    lat1->SetTextColor(1);
    lat1->DrawLatex(0.17, 0.2, "4.2 < Q^{2} [GeV] < 5");
    c1->Print(Form("Figs/Qp2_Scan2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Scan2_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Scan2_Run_%d.root", run));

    h_xi_xxGPD2->Draw();
    for (int i = 0; i < nQp2Bins; i++) {
        h_xi_xxGPD5_[i] = (TH2D*) file_in->Get(Form("h_xi_xxGPD5_%d", i));
        h_xi_xxGPD5_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD5_[i]->Draw("Same");
        //leg1->AddEntry(h_xi_xxGPD2_[i], Form("%1.2f < Q^{2}[GeV^{2}] < %1.2f", i+1., i+2.) );
        lat1->SetTextColor(i + 2);
        lat1->DrawLatex(0.72, 0.12 + i * 0.03, Form("%1.2f < Q^{'2}[GeV^{2}] < %1.2f", i + 1., i + 2.));
    }
    lat1->SetTextColor(1);
    lat1->DrawLatex(0.17, 0.2, "5. < Q^{2} [GeV] < 6.");
    c1->Print(Form("Figs/Qp2_Scan3_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Scan3_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Scan3_Run_%d.root", run));

    h_xi_xxGPD2->Draw();
    for (int i = 0; i < nQp2Bins; i++) {
        h_xi_xxGPD3_[i] = (TH2D*) file_in->Get(Form("h_xi_xxGPD3_%d", i));
        h_xi_xxGPD3_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD3_[i]->Draw("Same");
        //leg1->AddEntry(h_xi_xxGPD2_[i], Form("%1.2f < Q^{2}[GeV^{2}] < %1.2f", i+1., i+2.) );
        lat1->SetTextColor(i + 2);
        if (i >= 2 && i <= 5) {
            lat1->DrawLatex(0.72, 0.12 + i * 0.03, Form("%1.2f < Q^{2}[GeV^{2}] < %1.2f", i + 1., i + 2.));
        }
    }
    lat1->SetTextColor(1);
    lat1->DrawLatex(0.17, 0.2, "3. < Q^{'2} [GeV] < 4.");
    c1->Print(Form("Figs/Qp2_Scan4_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Scan4_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Scan5_Run_%d.root", run));

    c1->SetTopMargin(0.07);
    const double xB_max1 = 0.3;
    const double xB_min1 = 0.2;

    TFile *file_in12GeV = new TFile(Form("AnaGrape_%d.root", run_11GeV));
    c1->SetLogz();
    TH2D *h_Q2_xB2_11GeV = (TH2D*) file_in12GeV->Get("h_Q2_xB2");
    h_Q2_xB2_11GeV->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2_11GeV->Draw("colz");
    f_Q2_FixW->SetFillStyle(3004);
    f_Q2_FixW->SetFillColor(1);
    f_Q2_FixW->Draw("Same");
    c1->Print(Form("Figs/Q2_xB_11GeV_Run_%d.pdf", run_11GeV));
    c1->Print(Form("Figs/Q2_xB_11GeV_Run_%d.png", run_11GeV));
    c1->Print(Form("Figs/Q2_xB_11GeV_Run_%d.root", run_11GeV));

    int binxB_max1 = h_Q2_xB2_11GeV->GetXaxis()->FindBin(xB_max1);
    int binxB_min1 = h_Q2_xB2_11GeV->GetXaxis()->FindBin(xB_min1);
    TH1D *h_Q2Proj_11GeV = (TH1D*) h_Q2_xB2_11GeV->ProjectionY("h_Q2Proj_11GeV", binxB_min1, binxB_max1);
    h_Q2Proj_11GeV->Scale(1. / h_Q2Proj_11GeV->GetBinWidth(10));
    h_Q2Proj_11GeV->SetTitle("; Q^{2} [GeV^{2}]; dN/dQ^{2} [1/GeV^{2}]");
    h_Q2Proj_11GeV->Draw();



    TH2D *h_Q2_xB2_22GeV = (TH2D*) file_in->Get("h_Q2_xB2");
    h_Q2_xB2_22GeV->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2_22GeV->Draw("colz");
    c1->Print(Form("Figs/Q2_xB_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Q2_xB_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Q2_xB_22GeV_Run_%d.root", run));

    binxB_max1 = h_Q2_xB2_22GeV->GetXaxis()->FindBin(xB_max1);
    binxB_min1 = h_Q2_xB2_22GeV->GetXaxis()->FindBin(xB_min1);
    TH1D *h_Q2Proj_22GeV = (TH1D*) h_Q2_xB2_22GeV->ProjectionY("h_Q2Proj_22GeV", binxB_min1, binxB_max1);
    h_Q2Proj_22GeV->Scale(1. / h_Q2Proj_22GeV->GetBinWidth(10));
    h_Q2Proj_22GeV->SetTitle("; Q^{2} [GeV^{2}]; dN/dQ^{2} [1/GeV^{2}]");
    h_Q2Proj_22GeV->SetLineColor(2);
    h_Q2Proj_22GeV->Draw();

    h_Q2Proj_11GeV->Draw("");
    h_Q2Proj_22GeV->Draw("Same");
    lat1->SetTextSize(0.04);
    lat1->SetTextColor(4);
    lat1->DrawLatex(0.5, 0.8, "11 GeV beam");
    lat1->SetTextColor(2);
    lat1->DrawLatex(0.5, 0.7, "22 GeV beam");
    c1->Print("Figs/Q2_Dep_Compare_22Vs11GeV.pdf");
    c1->Print("Figs/Q2_Dep_Compare_22Vs11GeV.png");
    c1->Print("Figs/Q2_Dep_Compare_22Vs11GeV.root");


    //=======
    c1->SetLeftMargin(0.15);
    TH2D *h_Q2_xB2_11GeV_2 = (TH2D*) file_in12GeV->Get("h_Q2_xB3");

    h_Q2_xB2_11GeV_2->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2_11GeV_2->Draw("colz");
    c1->Print(Form("Figs/Q2_xB_11GeV_2_Run_%d.pdf", run_11GeV));
    c1->Print(Form("Figs/Q2_xB_11GeV_2_Run_%d.png", run_11GeV));
    c1->Print(Form("Figs/Q2_xB_11GeV_2_Run_%d.root", run_11GeV));

    binxB_max1 = h_Q2_xB2_11GeV_2->GetXaxis()->FindBin(xB_max1);
    binxB_min1 = h_Q2_xB2_11GeV_2->GetXaxis()->FindBin(xB_min1);
    TH1D *h_Q2Proj_11GeV_2 = (TH1D*) h_Q2_xB2_11GeV_2->ProjectionY("h_Q2Proj_11GeV_2", binxB_min1, binxB_max1);
    h_Q2Proj_11GeV_2->Scale(1. / h_Q2Proj_11GeV->GetBinWidth(10));
    h_Q2Proj_11GeV_2->SetTitle("; Q^{2} [GeV^{2}]; dN/dQ^{2} [1/GeV^{2}]");
    h_Q2Proj_11GeV_2->Draw();



    TH2D *h_Q2_xB2_22GeV_2 = (TH2D*) file_in->Get("h_Q2_xB3");
    h_Q2_xB2_22GeV_2->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2_22GeV_2->Draw("colz");
    c1->Print(Form("Figs/Q2_xB_22GeV_2_Run_%d.pdf", run));
    c1->Print(Form("Figs/Q2_xB_22GeV_2_Run_%d.png", run));
    c1->Print(Form("Figs/Q2_xB_22GeV_2_Run_%d.root", run));

    binxB_max1 = h_Q2_xB2_22GeV_2->GetXaxis()->FindBin(xB_max1);
    binxB_min1 = h_Q2_xB2_22GeV_2->GetXaxis()->FindBin(xB_min1);
    TH1D *h_Q2Proj_22GeV_2 = (TH1D*) h_Q2_xB2_22GeV_2->ProjectionY("h_Q2Proj_22GeV_2", binxB_min1, binxB_max1);
    h_Q2Proj_22GeV_2->Scale(1. / h_Q2Proj_22GeV_2->GetBinWidth(10));
    h_Q2Proj_22GeV_2->SetTitle("; Q^{2} [GeV^{2}]; dN/dQ^{2} [1/GeV^{2}]");
    h_Q2Proj_22GeV_2->SetLineColor(2);
    h_Q2Proj_22GeV->Draw();

    h_Q2Proj_11GeV_2->Draw("");
    h_Q2Proj_22GeV_2->Draw("Same");
    lat1->SetTextSize(0.04);
    lat1->SetTextColor(4);
    lat1->DrawLatex(0.5, 0.8, "11 GeV beam");
    lat1->SetTextColor(2);
    lat1->DrawLatex(0.5, 0.7, "22 GeV beam");
    c1->Print("Figs/Q2_Dep_Compare_22Vs11GeV_2.pdf");
    c1->Print("Figs/Q2_Dep_Compare_22Vs11GeV_2.png");
    c1->Print("Figs/Q2_Dep_Compare_22Vs11GeV_2.root");


    // ******* New kinematic bins, that we chose with Stepan in July 2024 for the 22 GeV Discusson meeting.
    c1->SetLogz(0);
    const double pol = 0.8;
    
    //TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);
    lat1->SetTextColor(4);
    lat1->SetTextSize(0.04);

    TLine *line1 = new TLine();
    line1->SetLineColor(95);
    line1->SetLineWidth(3);
    
    const int Qp2_Max = 5.;
    const int Qp2_Min = 4.;
    const int Q2_Max = 5;
    const int Q2_Min = 4;

    const int nQ2bins_new = 4;
    const int nQp2bins_new = 4;
    const double Q2_edges_new[nQ2bins_new + 1] = {2., 3.5, 5., 7., 10.};
    const double Qp2_edges_new[nQp2bins_new + 1] = {1.44, 3., 4., 6., 10.};
    
    const double tMcut = 0.4;

    h_xB_tM2->Draw("colz");
    line1->DrawLine(tMcut, 0., tMcut, 0.6);
    c1->Print(Form("Figs/xB_tM2new_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/xB_tM2new_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/xB_tM2new_22GeV_Run_%d.root", run));

    //TH2D *h_Qp2_vs_Q2_6 = (TH2D*) file_in->Get("h_Qp2_vs_Q2_6");
    h_Qp2_vs_Q2_6->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2}[GeV^{2}]");
    h_Qp2_vs_Q2_6->Draw("colz");

    line1->DrawLine(Q2_edges_new[0], Qp2_Max, Q2_edges_new[nQ2bins_new], Qp2_Max);
    line1->DrawLine(Q2_edges_new[0], Qp2_Min, Q2_edges_new[nQ2bins_new], Qp2_Min);
    for (int i = 0; i < nQ2bins_new + 1; i++) {
        line1->DrawLine(Q2_edges_new[i], Qp2_Min, Q2_edges_new[i], Qp2_Max);
    }

    c1->Print(Form("Figs/Qp2_vs_Q2_6_Q2Scan_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Q2Scan_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Q2Scan_22GeV_Run_%d.root", run));

    h_Qp2_vs_Q2_6->Draw("colz");
    line1->SetLineColor(6);
    line1->SetLineWidth(3);
    line1->DrawLine(Q2_Min, Qp2_edges_new[0], Q2_Min, Qp2_edges_new[nQp2bins_new]);
    line1->DrawLine(Q2_Max, Qp2_edges_new[0], Q2_Max, Qp2_edges_new[nQp2bins_new]);
    for (int i = 0; i < nQ2bins_new + 1; i++) {
        line1->DrawLine(Q2_Min, Qp2_edges_new[i], Q2_Max, Qp2_edges_new[i]);
    }
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Qp2Scan_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Qp2Scan_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_vs_Q2_6_Qp2Scan_22GeV_Run_%d.root", run));


    TH2D * h_xi_xxGPD_newQ2Scan_[nQ2bins_new];
    TH1D * h_Phi_LH_newQ2Scan_[nQ2bins_new];
    TH2D * h_Qp2_Q2_newQ2Scan_[nQ2bins_new];
    TH2D * h_tM_xB_newQ2Scan_[nQ2bins_new];
    TGraphErrors * gr_A_LH_newQ2Scan_[nQ2bins_new];

    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < nQ2bins_new; i++) {
        h_xi_xxGPD_newQ2Scan_[i] = (TH2D*) file_in->Get(Form("h_xi_xxGPD_newQ2Scan_%d", i));
        h_xi_xxGPD_newQ2Scan_[i]->SetMarkerColor(i + 2);
        lat1->SetTextColor(i + 2);
        h_xi_xxGPD_newQ2Scan_[i]->Draw("Same scat");
        lat1->DrawLatex(0.65, 0.4 - 0.05 * i, Form("Q^{2}#in (%1.1f-%1.1f)GeV", Q2_edges_new[i], Q2_edges_new[i + 1]));
    }
    c1->Print(Form("Figs/xi_vs_x_New_Q2Scan_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_New_Q2Scan_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_New_Q2Scan_22GeV_Run_%d.root", run));

    ofstream out_Q2Scan(Form("Dumps/Q2_Scan_Kinematics_22GeV_%d.dat", run));
    out_Q2Scan << "bin" << setw(10) << "Q2" << setw(12) << "Qp2" << setw(12) << "-t" << setw(12) << "xB" << endl;

    lat1->SetTextColor(4);
    for (int i = 0; i < nQ2bins_new; i++) {
        h_Phi_LH_newQ2Scan_[i] = (TH1D*) file_in->Get(Form("h_Phi_LH_newQ2Scan_%d", i));
        h_Phi_LH_newQ2Scan_[i]->SetMinimum(0);
        h_Phi_LH_newQ2Scan_[i]->SetTitle("; #phi_{LH} [deg]");
        h_Phi_LH_newQ2Scan_[i]->Draw();
        c1->Print(Form("Figs/Phi_LH_Q2Scan_Bin_%d_22GeV_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Phi_LH_Q2Scan_Bin_%d_22GeV_Run_%d.png", i, run));
        c1->Print(Form("Figs/Phi_LH_Q2Scan_Bin_%d_22GeV_Run_%d.root", i, run));
        ofstream histDump(Form("Dumps/Phi_LH_Q2Scan_Bin_%d_22GeV_Run_%d.dat", i, run));
        ifstream inp_asym(Form("Dumps/Asym_Q2Scan_Bin%d_22GeV.dat", i));

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

        h_Qp2_Q2_newQ2Scan_[i] = (TH2D*) file_in->Get(Form("h_Qp2_Q2_newQ2Scan_%d", i));
        double Q2 = h_Qp2_Q2_newQ2Scan_[i]->GetMean(1);
        double Qp2 = h_Qp2_Q2_newQ2Scan_[i]->GetMean(2);
        h_tM_xB_newQ2Scan_[i] = (TH2D*) file_in->Get(Form("h_tM_xB_newQ2Scan_%d", i));
        double tM = h_tM_xB_newQ2Scan_[i]->GetMean(1);
        double xB = h_tM_xB_newQ2Scan_[i]->GetMean(2);
        out_Q2Scan << i << setw(12) << Q2 << setw(12) << Qp2 << setw(12) << tM << setw(12) << xB << endl;

        gr_A_LH_newQ2Scan_[i]->Draw("AP");
        lat1->DrawLatex(0.04, 0.93, Form("Q^{2}=%1.2f GeV^{2}  Q^{'2} = %1.2f GeV^{2}  xB=%1.2f  -t=%1.2f GeV^{2}", Q2, Qp2, xB, tM));

        c1->Print(Form("Figs/Asym_Q2Scan_Bin_%d_22GeV_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Asym_Q2Scan_Bin_%d_22GeV_Run_%d.png", i, run));
        c1->Print(Form("Figs/Asym_Q2Scan_Bin_%d_22GeV_Run_%d.root", i, run));
    }

    TH2D * h_xi_xxGPD_newQp2Scan_[nQ2bins_new];
    TH1D * h_Phi_LH_newQp2Scan_[nQ2bins_new];
    TH2D * h_Qp2_Q2_newQp2Scan_[nQ2bins_new];
    TH2D * h_tM_xB_newQp2Scan_[nQ2bins_new];
    TGraphErrors * gr_A_LH_newQp2Scan_[nQ2bins_new];

    ofstream out_Qp2Scan(Form("Dumps/Qp2_Scan_Kinematics_22GeV_%d.dat", run));

    h_xi_xxGPD2->Draw("Scat");
    for (int i = 0; i < nQp2bins_new; i++) {
        h_xi_xxGPD_newQp2Scan_[i] = (TH2D*) file_in->Get(Form("h_xi_xxGPD_newQp2Scan_%d", i));
        h_xi_xxGPD_newQp2Scan_[i]->SetMarkerColor(i + 2);
        h_xi_xxGPD_newQp2Scan_[i]->Draw("Same scat");
        lat1->SetTextColor(i + 2);
        lat1->DrawLatex(0.65, 0.4 - 0.05 * i, Form("Q^{'2}#in (%1.1f-%1.1f)GeV", Qp2_edges_new[i], Qp2_edges_new[i + 1]));
    }
    c1->Print(Form("Figs/xi_vs_x_New_Qp2Scan_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/xi_vs_x_New_Qp2Scan_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/xi_vs_x_New_Qp2Scan_22GeV_Run_%d.root", run));

    lat1->SetTextColor(4);
    out_Qp2Scan << "bin" << setw(10) << "Q2" << setw(12) << "Qp2" << setw(12) << "-t" << setw(12) << "xB" << endl;
    for (int i = 0; i < nQp2bins_new; i++) {
        h_Phi_LH_newQp2Scan_[i] = (TH1D*) file_in->Get(Form("h_Phi_LH_newQp2Scan_%d", i));
        h_Phi_LH_newQp2Scan_[i]->SetMinimum(0);
        h_Phi_LH_newQp2Scan_[i]->SetTitle("; #phi_{LH} [deg]");
        h_Phi_LH_newQp2Scan_[i]->Draw();
        c1->Print(Form("Figs/Phi_LH_Qp2Scan_Bin_%d_22GeV_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Phi_LH_Qp2Scan_Bin_%d_22GeV_Run_%d.png", i, run));
        c1->Print(Form("Figs/Phi_LH_Qp2Scan_Bin_%d_22GeV_Run_%d.root", i, run));

        ofstream histDump(Form("Dumps/Phi_LH_Qp2Scan_Bin_%d_22GeV_Run_%d.dat", i, run));
        ifstream inp_asym(Form("Dumps/Asym_Qp2Scan_Bin%d_22GeV.dat", i));

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

        h_Qp2_Q2_newQp2Scan_[i] = (TH2D*) file_in->Get(Form("h_Qp2_Q2_newQp2Scan_%d", i));
        double Q2 = h_Qp2_Q2_newQp2Scan_[i]->GetMean(1);
        double Qp2 = h_Qp2_Q2_newQp2Scan_[i]->GetMean(2);
        h_tM_xB_newQp2Scan_[i] = (TH2D*) file_in->Get(Form("h_tM_xB_newQp2Scan_%d", i));
        double tM = h_tM_xB_newQp2Scan_[i]->GetMean(1);
        double xB = h_tM_xB_newQp2Scan_[i]->GetMean(2);
        out_Qp2Scan << i << setw(12) << Q2 << setw(12) << Qp2 << setw(12) << tM << setw(12) << xB << endl;

        gr_A_LH_newQ2Scan_[i]->Draw("AP");
        lat1->DrawLatex(0.04, 0.93, Form("Q^{2}=%1.2f GeV^{2}  Q^{'2} = %1.2f GeV^{2}  xB=%1.2f  -t=%1.2f GeV^{2}", Q2, Qp2, xB, tM));

        c1->Print(Form("Figs/Asym_Qp2Scan_Bin_%d_22GeV_Run_%d.pdf", i, run));
        c1->Print(Form("Figs/Asym_Qp2Scan_Bin_%d_22GeV_Run_%d.png", i, run));
        c1->Print(Form("Figs/Asym_Qp2Scan_Bin_%d_22GeV_Run_%d.root", i, run));

    }

    return 0;
}