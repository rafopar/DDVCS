/* 
 * File:   DrawPlots.cc
 * Author: rafopar
 *
 * Created on October 2, 2024, 8:13 PM
 */

#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

template <typename T> void FormatHist(T h);

/*
 * 
 */
int main(int argc, char** argv) {

    if( argc < 2 ){
        cout<<"Please provide the run number"<<endl;
        cout<<"Exiting"<<endl;
        exit(1);
    }
    
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    
    int run = atoi(argv[1]);
    
    TCanvas c1("c1", "", 950, 950);
    c1.SetTopMargin(0.04);
    c1.SetRightMargin(0.04);
    c1.SetLeftMargin(0.13);
    c1.SetBottomMargin(0.13);
    
    TLatex lat1;
    lat1.SetNDC();
    
    TFile file_in(Form("AnaJPsi_Run_%d.root", run), "Read");
    
    
    TH2D *h_th_P_mup1 = (TH2D*)file_in.Get("h_th_P_mup1");
    h_th_P_mup1->SetTitle("; P_{#mu^{+}} [GeV]; #theta_{#mu^{+}} ");
    FormatHist(h_th_P_mup1);
    h_th_P_mup1->Draw();
    c1.Print(Form("Figs/th_P_mup1_Run_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mup1_Run_%d.png", run));
    c1.Print(Form("Figs/th_P_mup1_Run_%d.root", run));

    TH2D *h_th_P_mum1 = (TH2D*)file_in.Get("h_th_P_mum1");
    h_th_P_mum1->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) ");
    FormatHist(h_th_P_mum1);
    h_th_P_mum1->Draw();
    c1.Print(Form("Figs/th_P_mum1_Run_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mum1_Run_%d.png", run));
    c1.Print(Form("Figs/th_P_mum1_Run_%d.root", run));
        
    TH1D *h_Mmumu_Rec1 = (TH1D*)file_in.Get("h_Mmumu_Rec1");
    h_Mmumu_Rec1->SetLineColor(1);
    h_Mmumu_Rec1->SetLineWidth(2);    
    h_Mmumu_Rec1->SetTitle("; M(#mu#mu) [GeV]");
    FormatHist(h_Mmumu_Rec1);
    
    TH1D *h_Mmumu_Rec_Corr1 = (TH1D*)file_in.Get("h_Mmumu_Rec_Corr1");
    h_Mmumu_Rec_Corr1->SetTitle("; M(#mu#mu) [GeV]");
    h_Mmumu_Rec_Corr1->SetLineColor(2);
    h_Mmumu_Rec_Corr1->SetLineWidth(2);
    FormatHist(h_Mmumu_Rec_Corr1);
    
    TH1D *h_Mmumu_Rec_Corr2 = (TH1D*)file_in.Get("h_Mmumu_Rec_Corr2");
    h_Mmumu_Rec_Corr2->SetTitle("; M(#mu#mu) [GeV]");
    h_Mmumu_Rec_Corr2->SetLineColor(6);
    h_Mmumu_Rec_Corr2->SetLineWidth(2);
    FormatHist(h_Mmumu_Rec_Corr2);
    
    TH1D *h_Mmumu_Rec_Corr3 = (TH1D*)file_in.Get("h_Mmumu_Rec_Corr3");
    h_Mmumu_Rec_Corr3->SetTitle("; M(#mu#mu) [GeV]");
    h_Mmumu_Rec_Corr3->SetLineColor(4);
    h_Mmumu_Rec_Corr3->SetLineWidth(2);
    FormatHist(h_Mmumu_Rec_Corr3);
    
    TH1D *h_Mmumu_Rec_AngleFixCorr1 = (TH1D*)file_in.Get("h_Mmumu_Rec_AngleFixCorr1");
    h_Mmumu_Rec_AngleFixCorr1->SetTitle("; M(#mu#mu) [GeV]");
    h_Mmumu_Rec_AngleFixCorr1->SetLineColor(8);
    h_Mmumu_Rec_AngleFixCorr1->SetLineWidth(2);
    FormatHist(h_Mmumu_Rec_AngleFixCorr1);
    
    c1.SetGridx();
    TLegend *leg1 = new TLegend(0.15, 0.7, 0.5, 0.95);
    leg1->SetBorderSize(0);
    leg1->AddEntry(h_Mmumu_Rec1, "Uncorrected");
    leg1->AddEntry(h_Mmumu_Rec_Corr1, "Mom. corrected");
    leg1->AddEntry(h_Mmumu_Rec_Corr2, "Mom & #theta corrected");
    leg1->AddEntry(h_Mmumu_Rec_Corr3, "Mom & #theta & #phi corrected");
    leg1->AddEntry(h_Mmumu_Rec_AngleFixCorr1, "Mom corrected  #theta & #phi Fixed");
    h_Mmumu_Rec1->Draw();
    h_Mmumu_Rec_Corr1->Draw("Same");
    double max = TMath::Max( h_Mmumu_Rec1->GetMaximum(), h_Mmumu_Rec_Corr1->GetMaximum() );
    h_Mmumu_Rec_Corr2->Draw("Same");
    max = TMath::Max(max, h_Mmumu_Rec_Corr2->GetMaximum() );
    h_Mmumu_Rec_Corr3->Draw("Same");
    max = TMath::Max(max, h_Mmumu_Rec_Corr3->GetMaximum() );
    h_Mmumu_Rec_AngleFixCorr1->Draw("Same");
    max = TMath::Max(max, h_Mmumu_Rec_AngleFixCorr1->GetMaximum() );
    h_Mmumu_Rec1->SetMaximum(1.05*max);
    
    leg1->Draw();
    double rms_rec = h_Mmumu_Rec1->GetRMS();
    double rms_cor1 = h_Mmumu_Rec_Corr1->GetRMS();
    double rms_cor2 = h_Mmumu_Rec_Corr2->GetRMS();
    double rms_cor3 = h_Mmumu_Rec_Corr3->GetRMS();
    double rms_AngleFixCor1 = h_Mmumu_Rec_AngleFixCorr1->GetRMS();
    lat1.SetTextSize(0.03);
    lat1.SetTextColor(1);
    lat1.DrawLatex(0.15, 0.6, Form("Uncor, RMS = %1.0f MeV", rms_rec*1000));
    lat1.SetTextColor(2);
    lat1.DrawLatex(0.15, 0.55, Form("p cor, RMS = %1.0f MeV", rms_cor1*1000));
    lat1.SetTextColor(6);
    lat1.DrawLatex(0.15, 0.5, Form("p + #theta cor, RMS = %1.0f MeV", rms_cor2*1000));
    lat1.SetTextColor(4);
    lat1.DrawLatex(0.15, 0.45, Form("p + #theta + #phi cor, RMS = %1.0f MeV", rms_cor3*1000));
    lat1.SetTextColor(8);
    lat1.DrawLatex(0.15, 0.4, Form("p corr + #theta + #phi Fix, RMS = %1.0f MeV", rms_AngleFixCor1*1000));
    c1.Print(Form("Figs/MmumuCorretions_Run_%d.pdf", run));
    c1.Print(Form("Figs/MmumuCorretions_Run_%d.png", run));
    c1.Print(Form("Figs/MmumuCorretions_Run_%d.root", run));
            
    return 0;
}

template <typename T>
void FormatHist(T h){

    //TH2D *h = new TH2D("dsd", "", 100, 0, 100, 200, 0, 100);
    
    h->SetTitleSize(0.05, "X");
    h->SetTitleSize(0.05, "Y");
    h->SetLabelSize(0.05, "X");
    h->SetLabelSize(0.05, "Y");
}