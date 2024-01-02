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
int DrawPlots_22GeV() {

    gStyle->SetOptStat(0);

    const int run = 4;
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

    TH2D *h_Qp2_vs_Q2_3 = (TH2D*) file_in->Get("h_Qp2_vs_Q2_3");
    h_Qp2_vs_Q2_3->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}] ");
    h_Qp2_vs_Q2_3->SetAxisRange(2.5, 8.5, "X");
    h_Qp2_vs_Q2_3->Draw("colz");
    c1->Print(Form("Figs/Qp2_Q2_3_22GeV_Run_%d.pdf", run));
    c1->Print(Form("Figs/Qp2_Q2_3_22GeV_Run_%d.png", run));
    c1->Print(Form("Figs/Qp2_Q2_3_22GeV_Run_%d.root", run));



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


    return 0;
}

