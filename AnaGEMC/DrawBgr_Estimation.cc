/* 
 * File:   DrawBgr_Estimation.cc
 * Author: rafopar
 *
 * Created on March 21, 2025, 3:30 PM
 */

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <cstdlib>

using namespace std;

template <typename T> void FormatHist(T h);

/*
 * 
 */
void DrawBgr_Estimation() {

    gStyle->SetOptStat(0);
    
    /*
     *  After hadd, we need to scale rates by 1./nRuns, so that is why we need nRuns
     */
    const int nRuns = 57;

    const int nDays = 200;
    const int secPerDay = 3600 * 24;
    const int TotTimeinSec = nDays*secPerDay;

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);
    
    //TFile *file_pionBgr = new TFile("AnaDDVCS_Run_RGA_8Runs.root", "Read");
    //TFile *file_pionBgr = new TFile("AnaDDVCS_Run_RGA_7Runs_Lead.root", "Read");
    //TFile *file_pionBgr = new TFile("AnaDDVCS_Run_RGA_9Runs_Lead_55cm.root", "Read");
//    TFile *file_pionBgr = new TFile("AnaDDVCS_Run_RGA_10Runs_Lead_60cm.root", "Read");
    TFile *file_pionBgr = new TFile("AnaDDVCS_Run_RGA_57Runs_Lead_60cm.root", "Read");

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetLogy();
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    c1->SetLeftMargin(0.18);

    const int nStages = 5;
    /* 
     * 0: All events
     * 1: Has at least muon pair detected
     * 2: Has at least 1mu-, 1mu+ and an electron
     * 3: Pass the the missing mass cut
     * 4: Pass the -t cut
     * 5: Both muons have P > 2 GeV
     */

    int cols_[nStages] = {1, 8, 2, 4, 6};
    int ind_[nStages]= {1, 6, 2, 3, 4};
    
    TH1D * h_Mmumu_MC[nStages];

    for (int i = 0; i < nStages; i++) {
        //h_Mmumu_MC[i] = (TH1D*) file_pionBgr->Get(Form("h_Mmumu_MC%d", i+1));
        h_Mmumu_MC[i] = (TH1D*) file_pionBgr->Get(Form("h_Mmumu_MC%d", ind_[i]));
        h_Mmumu_MC[i]->Scale(1. / double(nRuns));
        h_Mmumu_MC[i]->Scale(1. / double(TotTimeinSec));
        double binWidth = h_Mmumu_MC[i]->GetBinWidth(10);
        h_Mmumu_MC[i]->SetTitle(Form("; Invariant mass [GeV]; Rate [Hz/%1.2f GeV]", binWidth));
        
        double tot_Rate = h_Mmumu_MC[i]->Integral();
        
        FormatHist(h_Mmumu_MC[i]);
        h_Mmumu_MC[i]->SetLineColor(cols_[i]);
        if (i == 0) {
            h_Mmumu_MC[i]->SetMinimum(1.e-4);
            h_Mmumu_MC[i]->Draw("hist");
        } else {
            h_Mmumu_MC[i]->Draw("hist Same");
        }
        
        lat1->SetTextColor(cols_[i]);
        lat1->DrawLatex(0.25, 0.6 - 0.05*i, Form("Tot. Rate = %1.4f Hz", tot_Rate));
    }
    
    c1->Print("Figs/pion_pair_background.pdf");
    c1->Print("Figs/pion_pair_background.png");
    c1->Print("Figs/pion_pair_background.root");
}

template <typename T>
void FormatHist(T h) {

    h->SetLineWidth(3);
    h->SetTitleSize(0.05, "X");
    h->SetTitleSize(0.05, "Y");
    h->SetLabelSize(0.05, "X");
    h->SetLabelSize(0.05, "Y");
}
