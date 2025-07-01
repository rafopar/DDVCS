/* 
 * File:   Draw_ISR_Effects.cc
 * Author: rafopar
 *
 * Created on May 1, 2025, 9:38 AM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
void Draw_ISR_Effects() {

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);

    TFile *file_ISR = new TFile("AnaDDVCS_Run_32005.root", "Read");
    TFile *file_NOISR = new TFile("AnaDDVCS_Run_33005.root", "Read");

    TH2D *h_Qp2_Q2_xi_x_0_1_ISR = (TH2D*) file_ISR->Get("h_Qp2_Q2_xi_x_0_1");
    TH1D *h_Q2_xi_x_0_1_ISR = (TH1D*) h_Qp2_Q2_xi_x_0_1_ISR->ProjectionX("h_Q2_xi_x_0_1_ISR", 1, h_Qp2_Q2_xi_x_0_1_ISR->GetNbinsY());
    h_Q2_xi_x_0_1_ISR->SetTitle("; Q^{2} [GeV^{2}");
    h_Q2_xi_x_0_1_ISR->SetLineColor(4);
    TH1D *h_Q2_xi_x_0_1_ISR_RB = (TH1D*) h_Q2_xi_x_0_1_ISR->Rebin(5, "h_Q2_xi_x_0_1_ISR_RB");

    TH2D *h_Qp2_Q2_xi_x_0_1_NOISR = (TH2D*) file_NOISR->Get("h_Qp2_Q2_xi_x_0_1");
    TH1D *h_Q2_xi_x_0_1_NOISR = (TH1D*) h_Qp2_Q2_xi_x_0_1_NOISR->ProjectionX("h_Q2_xi_x_0_1_NOISR", 1, h_Qp2_Q2_xi_x_0_1_NOISR->GetNbinsY());
    h_Q2_xi_x_0_1_NOISR->SetTitle("; Q^{2} [GeV^{2}");
    h_Q2_xi_x_0_1_NOISR->SetLineColor(2);
    TH1D *h_Q2_xi_x_0_1_NOISR_RB = (TH1D*) h_Q2_xi_x_0_1_NOISR->Rebin(5, "h_Q2_xi_x_0_1_NOISR_RB");


    h_Q2_xi_x_0_1_NOISR_RB->Draw();
    h_Q2_xi_x_0_1_ISR_RB->Draw("Same");

    TRatioPlot *rp1 = new TRatioPlot(h_Q2_xi_x_0_1_ISR_RB, h_Q2_xi_x_0_1_NOISR_RB);
    rp1->Draw();

    TH1D *h_Q2_1_ISR = (TH1D*) file_ISR->Get("h_Q2_1");
    h_Q2_1_ISR->SetTitle("; Q^{2} [GeV^{2}]");
    h_Q2_1_ISR->SetLineColor(4);

    TH1D *h_Q2_1_NOISR = (TH1D*) file_NOISR->Get("h_Q2_1");
    h_Q2_1_NOISR->SetTitle("; Q^{2} [GeV^{2}]");
    h_Q2_1_NOISR->SetLineColor(2);

    TRatioPlot *rp2 = new TRatioPlot(h_Q2_1_ISR, h_Q2_1_NOISR);
    rp2->Draw();
    rp2->GetLowerRefGraph()->SetMaximum(1.3);
    rp2->GetLowerRefGraph()->SetMinimum(0.9);
    c1->Print("Figs/ISR_Effect_Q2_1.pdf");
    c1->Print("Figs/ISR_Effect_Q2_1.png");
    c1->Print("Figs/ISR_Effect_Q2_1.root");


    TH1D *h_Q2_MC_1_ISR = (TH1D*) file_ISR->Get("h_Q2_MC_1");
    h_Q2_MC_1_ISR->SetTitle("; Q^{2} [GeV^{2}]");
    h_Q2_MC_1_ISR->SetLineColor(4);

    TH1D *h_Q2_MC_1_NOISR = (TH1D*) file_NOISR->Get("h_Q2_MC_1");
    h_Q2_MC_1_NOISR->SetTitle("; Q^{2} [GeV^{2}]");
    h_Q2_MC_1_NOISR->SetLineColor(2);

    TRatioPlot *rp3 = new TRatioPlot(h_Q2_MC_1_ISR, h_Q2_MC_1_NOISR);
    rp3->Draw();
    rp3->GetLowerRefGraph()->SetMaximum(1.3);
    rp3->GetLowerRefGraph()->SetMinimum(0.9);
    c1->Print("Figs/ISR_Effect_Q2_MC_1.pdf");
    c1->Print("Figs/ISR_Effect_Q2_MC_1.png");
    c1->Print("Figs/ISR_Effect_Q2_MC_1.root");

}

