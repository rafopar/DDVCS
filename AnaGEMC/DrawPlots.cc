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
    
    TFile file_in(Form("AnaDDVCS_Run_%d.root", run), "Read");
    
    
    TH2D *h_th_P_mup1 = (TH2D*)file_in.Get("h_th_P_mup1");
    h_th_P_mup1->SetTitle("; P_{#mu^{+}} [GeV]; #theta_{#mu^{+}} ");
    FormatHist(h_th_P_mup1);
    h_th_P_mup1->Draw();
    c1.Print(Form("Figs/th_P_mup1_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mup1_%d.png", run));
    c1.Print(Form("Figs/th_P_mup1_%d.root", run));

    TH2D *h_th_P_mum1 = (TH2D*)file_in.Get("h_th_P_mum1");
    h_th_P_mum1->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-}) ");
    FormatHist(h_th_P_mum1);
    h_th_P_mum1->Draw();
    c1.Print(Form("Figs/th_P_mum1_%d.pdf", run));
    c1.Print(Form("Figs/th_P_mum1_%d.png", run));
    c1.Print(Form("Figs/th_P_mum1_%d.root", run));
    
    TH2D *h_MC_th_P_em2 = (TH2D*)file_in.Get("h_MC_th_P_em2");
    h_MC_th_P_em2->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-}) ");
    FormatHist(h_MC_th_P_em2);
    h_MC_th_P_em2->Draw();
    c1.Print(Form("Figs/th_P_em2_%d.pdf", run));
    c1.Print(Form("Figs/th_P_em2_%d.png", run));
    c1.Print(Form("Figs/th_P_em2_%d.root", run));
    
    TH1D *h_Mmis1 = (TH1D*)file_in.Get("h_Mmis1");
    h_Mmis1->SetLineColor(1);
    h_Mmis1->SetLineWidth(2);    
    h_Mmis1->SetTitle("; M_{X}^{2} [GeV^{2}]");
    FormatHist(h_Mmis1);
    
    TH1D *h_Mmis_corr1 = (TH1D*)file_in.Get("h_Mmis_corr1");
    h_Mmis_corr1->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_corr1->SetLineColor(2);
    h_Mmis_corr1->SetLineWidth(2);
    FormatHist(h_Mmis_corr1);
    
    TH1D *h_Mmis_corr2 = (TH1D*)file_in.Get("h_Mmis_corr2");
    h_Mmis_corr2->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_corr2->SetLineColor(6);
    h_Mmis_corr2->SetLineWidth(2);
    FormatHist(h_Mmis_corr2);
    
    TH1D *h_Mmis_corr3 = (TH1D*)file_in.Get("h_Mmis_corr3");
    h_Mmis_corr3->SetTitle("; M_{X}^{2} [GeV^{2}]");
    h_Mmis_corr3->SetLineColor(4);
    h_Mmis_corr3->SetLineWidth(2);
    FormatHist(h_Mmis_corr3);
    
    TLegend *leg1 = new TLegend(0.5, 0.65, 0.91, 0.9);
    leg1->SetBorderSize(0);
    leg1->AddEntry(h_Mmis1, "Uncorrected");
    leg1->AddEntry(h_Mmis_corr1, "Mom. corrected");
    leg1->AddEntry(h_Mmis_corr2, "Mom & #theta corrected");
    leg1->AddEntry(h_Mmis_corr3, "Mom & #theta & #phi corrected");
    h_Mmis_corr3->Draw();
    h_Mmis_corr2->Draw("Same");
    h_Mmis_corr1->Draw("Same");
    h_Mmis1->Draw("Same");
    leg1->Draw();
    c1.Print(Form("Figs/MmisCorretions_%d.pdf", run));
    c1.Print(Form("Figs/MmisCorretions_%d.png", run));
    c1.Print(Form("Figs/MmisCorretions_%d.root", run));
    
    TH1D *h_Delta_tM1 = (TH1D*)file_in.Get("h_Delta_tM1");
    h_Delta_tM1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM1->SetLineColor(4);
    h_Delta_tM1->SetLineWidth(2);
    FormatHist(h_Delta_tM1);
    
    TH1D *h_Delta_tM_Constrian1 = (TH1D*)file_in.Get("h_Delta_tM_Constrian1");
    h_Delta_tM_Constrian1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM_Constrian1->SetLineColor(2);
    h_Delta_tM_Constrian1->SetLineWidth(2);
    FormatHist(h_Delta_tM_Constrian1);
    
    TH1D *h_Delta_tM_FixMom_Constrain1 = (TH1D*)file_in.Get("h_Delta_tM_FixMom_Constrain1");
    h_Delta_tM_FixMom_Constrain1->SetTitle("; -(t - t_{MC}) [GeV^{2}]");
    h_Delta_tM_FixMom_Constrain1->SetLineColor(6);
    h_Delta_tM_FixMom_Constrain1->SetLineWidth(2);
    FormatHist(h_Delta_tM_FixMom_Constrain1);
        
    TH1D *h_Delta_tM_FixEnergy_Constrain1 = (TH1D*)file_in.Get("h_Delta_tM_FixEnergy_Constrain1");
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