/* 
 * File:   DrawPionRejection.cc
 * Author: rafopar
 *
 * Created on January 29, 2025, 6:15 PM
 */

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iostream>

#include <clas12AnaTools.h>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    gStyle->SetOptStat(0);

    if (argc < 2) {
        cout << "Please provide the Identifier as an argument" << endl;
        cout << "The program should run as ./AnaPionRejection.exe Identifier" << endl;
        cout << "As an example ./AnaPionRejection.exe piPlus_6GeV_22Deg_0Deg_5GeV_17deg_40deg" << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    std::string Identifier = argv[1];

    cout << "Identifier is " << Identifier << endl;

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.03);
    c1->SetRightMargin(0.13);
    c1->SetLeftMargin(0.12);

    TFile *file_in = new TFile(Form("pion_Rejection_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_MC1 = (TH2D*) file_in->Get("h_cosTh_P_pi_MC1");
    h_cosTh_P_pi_MC1->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_MC1->Draw();
    c1->Print(Form("Figs/cosTh_P_pi_MC1_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_MC1_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_MC1_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mup_MC2 = (TH2D*) file_in->Get("h_cosTh_P_pi_mup_MC2");
    h_cosTh_P_pi_mup_MC2->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_mup_MC2->Draw();
    c1->Print(Form("Figs/cosTh_P_pi_mup_MC2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_MC2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_MC2_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mup_Ratio = (TH2D*) h_cosTh_P_pi_mup_MC2->Clone("h_cosTh_P_pi_mup_Ratio");
    h_cosTh_P_pi_mup_Ratio->Divide(h_cosTh_P_pi_MC1);
    h_cosTh_P_pi_mup_Ratio->Scale(100.);
    h_cosTh_P_pi_mup_Ratio->SetMaximum(1.05);
    h_cosTh_P_pi_mup_Ratio->Draw("colz");
    c1->Print(Form("Figs/pion_mup_Survive_Efficiency_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mup_Survive_Efficiency_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mup_Survive_Efficiency_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mum_MC2 = (TH2D*) file_in->Get("h_cosTh_P_pi_mum_MC2");
    h_cosTh_P_pi_mum_MC2->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_mum_MC2->Draw();
    c1->Print(Form("Figs/cosTh_P_pi_mum_MC2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mum_MC2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mum_MC2_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mum_Ratio = (TH2D*) h_cosTh_P_pi_mum_MC2->Clone("h_cosTh_P_pi_mum_Ratio");
    h_cosTh_P_pi_mum_Ratio->Divide(h_cosTh_P_pi_MC1);
    h_cosTh_P_pi_mum_Ratio->Scale(100.);
    //h_cosTh_P_pi_mum_Ratio->SetMaximum(5.);
    h_cosTh_P_pi_mum_Ratio->Draw("colz");
    c1->Print(Form("Figs/pion_mum_Survive_Efficiency_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mum_Survive_Efficiency_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mum_Survive_Efficiency_%s.root", Identifier.c_str()));


    TH2D *h_cosTh_P_pi_mup_MC3 = (TH2D*) file_in->Get("h_cosTh_P_pi_mup_MC3");
    h_cosTh_P_pi_mup_MC3->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_mup_MC3->Draw();
    c1->Print(Form("Figs/cosTh_P_pi_mup_MC3_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_MC3_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_MC3_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mup_DirectSurvive = (TH2D*) h_cosTh_P_pi_mup_MC3->Clone("h_cosTh_P_pi_mup_DirectSurvive");
    h_cosTh_P_pi_mup_DirectSurvive->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_mup_DirectSurvive->Divide(h_cosTh_P_pi_MC1);
    h_cosTh_P_pi_mup_DirectSurvive->Scale(100.);
    h_cosTh_P_pi_mup_DirectSurvive->Draw("colz");
    c1->Print(Form("Figs/cosTh_P_pi_mup_DirectSurvive_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_DirectSurvive_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_DirectSurvive_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mup_InDirectSurvive = (TH2D*) h_cosTh_P_pi_mup_Ratio->Clone("h_cosTh_P_pi_mup_InDirectSurvive");
    h_cosTh_P_pi_mup_InDirectSurvive->Add(h_cosTh_P_pi_mup_DirectSurvive, -1.);
    h_cosTh_P_pi_mup_InDirectSurvive->SetMaximum(1);
    h_cosTh_P_pi_mup_InDirectSurvive->Draw("colz");
    c1->Print(Form("Figs/cosTh_P_pi_mup_InDirectSurvive_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_InDirectSurvive_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_InDirectSurvive_%s.root", Identifier.c_str()));

    
    TH2D *h_cosTh_P_pi_mup_Loose_MC2 = (TH2D*) file_in->Get("h_cosTh_P_pi_mup_Loose_MC2");
    h_cosTh_P_pi_mup_Loose_MC2->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_mup_Loose_MC2->Draw();
    c1->Print(Form("Figs/cosTh_P_pi_mup_Loose_MC2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_Loose_MC2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mup_Loose_MC2_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mup_Loose_Ratio = (TH2D*) h_cosTh_P_pi_mup_Loose_MC2->Clone("h_cosTh_P_pi_mup_Loose_Ratio");
    h_cosTh_P_pi_mup_Loose_Ratio->Divide(h_cosTh_P_pi_MC1);
    h_cosTh_P_pi_mup_Loose_Ratio->Scale(100.);
    h_cosTh_P_pi_mup_Loose_Ratio->SetMaximum(1.05);
    h_cosTh_P_pi_mup_Loose_Ratio->Draw("colz");
    c1->Print(Form("Figs/pion_mup_Loose_Survive_Efficiency_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mup_Loose_Survive_Efficiency_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mup_Loose_Survive_Efficiency_%s.root", Identifier.c_str()));
    
    TH2D *h_cosTh_P_pi_mum_Loose_MC2 = (TH2D*) file_in->Get("h_cosTh_P_pi_mum_Loose_MC2");
    h_cosTh_P_pi_mum_Loose_MC2->SetTitle("; P(MC) [GeV]; cos(#theta) (MC) ");
    h_cosTh_P_pi_mum_Loose_MC2->Draw();
    c1->Print(Form("Figs/cosTh_P_pi_mum_Loose_MC2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mum_Loose_MC2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_P_pi_mum_Loose_MC2_%s.root", Identifier.c_str()));

    TH2D *h_cosTh_P_pi_mum_Loose_Ratio = (TH2D*) h_cosTh_P_pi_mum_Loose_MC2->Clone("h_cosTh_P_pi_mum_Loose_Ratio");
    h_cosTh_P_pi_mum_Loose_Ratio->Divide(h_cosTh_P_pi_MC1);
    h_cosTh_P_pi_mum_Loose_Ratio->Scale(100.);
    h_cosTh_P_pi_mum_Loose_Ratio->SetMaximum(1.05);
    h_cosTh_P_pi_mum_Loose_Ratio->Draw("colz");
    c1->Print(Form("Figs/pion_mum_Loose_Survive_Efficiency_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mum_Loose_Survive_Efficiency_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pion_mum_Loose_Survive_Efficiency_%s.root", Identifier.c_str()));
    
    
    
    TH2D *h_cosTh_phi_pi_mup_MC2 = (TH2D*) file_in->Get("h_cosTh_phi_pi_mup_MC2");
    h_cosTh_phi_pi_mup_MC2->SetTitle("; #phi (MC) [deg]; cos(#theta) (MC) ");
    h_cosTh_phi_pi_mup_MC2->Draw();
    c1->Print(Form("Figs/cosTh_phi_pi_mup_MC2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_phi_pi_mup_MC2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/cosTh_phi_pi_mup_MC2_%s.root", Identifier.c_str()));

    

    TF1 *f_lin1 = new TF1("f_lin1", "[0]  + [1]*x", 12.);

    TH2D *h_PRec_PMC_mup1 = (TH2D*) file_in->Get("h_PRec_PMC_mup1");
    h_PRec_PMC_mup1->SetTitle("; P_{MC}(#pi^{+}) [GeV]; P_{Rec}(#pi^{+}) [GeV] ");
    h_PRec_PMC_mup1->Draw();
    f_lin1->SetParameters(-2.5, 1.);
    f_lin1->DrawCopy("Same");
    f_lin1->SetParameters(-0.8, 1.);
    f_lin1->DrawCopy("Same");
    c1->Print(Form("Figs/PRec_PMC_mup1_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/PRec_PMC_mup1_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/PRec_PMC_mup1_%s.root", Identifier.c_str()));

    TH2D *h_PRec_PMC_mum1 = (TH2D*) file_in->Get("h_PRec_PMC_mum1");
    h_PRec_PMC_mum1->SetTitle("; P_{MC}(#pi^{+}) [GeV]; P_{Rec}(#pi^{+}) [GeV] ");
    h_PRec_PMC_mum1->Draw();
    c1->Print(Form("Figs/PRec_PMC_mum1_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/PRec_PMC_mum1_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/PRec_PMC_mum1_%s.root", Identifier.c_str()));

    TF1 *f_SF_MIP = new TF1("f_SF_MIP", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 20.);
    f_SF_MIP->SetNpx(4500);

//    TH2D *h_SF_pos4 = (TH2D*)file_in->Get("h_SF_pos4");
//    h_SF_pos4->SetTitle("; P (positive) [GeV]; SF");
//    clas12AnaTools::SliceFit(h_SF_pos4, 0.5, 9., 30, f_SF_MIP, Form("SF_pos_%s", Identifier.c_str()));

    c1->cd();
    TH2D *h_SF_pos1 = (TH2D*)file_in->Get("h_SF_pos1");
    h_SF_pos1->SetTitle("; P (positive) [GeV]; SF");    
    h_SF_pos1->Draw();
    c1->Print(Form("Figs/SF_pos1_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/SF_pos1_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/SF_pos1_%s.root", Identifier.c_str()));
    
    TH2D *h_SF_pos5 = (TH2D*)file_in->Get("h_SF_pos5");
    h_SF_pos5->SetTitle("; P (positive) [GeV]; SF");    
    h_SF_pos5->Draw();
    c1->Print(Form("Figs/SF_pos5_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/SF_pos5_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/SF_pos5_%s.root", Identifier.c_str()));
    
    TH2D *h_SF_neg1 = (TH2D*)file_in->Get("h_SF_neg1");
    h_SF_neg1->SetTitle("; P (negative) [GeV]; SF");    
    h_SF_neg1->Draw();
    c1->Print(Form("Figs/SF_neg1_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/SF_neg1_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/SF_neg1_%s.root", Identifier.c_str()));
    
    TH2D *h_SF_neg5 = (TH2D*)file_in->Get("h_SF_neg5");
    h_SF_neg5->SetTitle("; P (negative) [GeV]; SF");    
    h_SF_neg5->Draw();
    c1->Print(Form("Figs/SF_neg5_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/SF_neg5_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/SF_neg5_%s.root", Identifier.c_str()));
    
    TH1D *h_pos_n_StripPCal2 = (TH1D*)file_in->Get("h_pos_n_StripPCal2");
    h_pos_n_StripPCal2->SetTitle("; Number of PCal Strips");
    h_pos_n_StripPCal2->Draw();
    c1->Print(Form("Figs/pos_N_PCalStrips2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pos_N_PCalStrips2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pos_N_PCalStrips2_%s.root", Identifier.c_str()));
    
    TH1D *h_pos_n_StripECin2 = (TH1D*)file_in->Get("h_pos_n_StripECin2");
    h_pos_n_StripECin2->SetTitle("; Number of EC_{in} Strips");
    h_pos_n_StripECin2->Draw();
    c1->Print(Form("Figs/pos_N_ECinStrips2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pos_N_ECinStrips2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pos_N_ECinStrips2_%s.root", Identifier.c_str()));
    
    TH1D *h_pos_n_StripECout2 = (TH1D*)file_in->Get("h_pos_n_StripECout2");
    h_pos_n_StripECout2->SetTitle("; Number of EC_{out} Strips");
    h_pos_n_StripECout2->Draw();
    c1->Print(Form("Figs/pos_N_ECoutStrips2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/pos_N_ECoutStrips2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/pos_N_ECoutStrips2_%s.root", Identifier.c_str()));
    
    TH1D *h_neg_n_StripPCal2 = (TH1D*)file_in->Get("h_neg_n_StripPCal2");
    h_neg_n_StripPCal2->SetTitle("; Number of PCal Strips");
    h_neg_n_StripPCal2->Draw();
    c1->Print(Form("Figs/neg_N_PCalStrips2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/neg_N_PCalStrips2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/neg_N_PCalStrips2_%s.root", Identifier.c_str()));
    
    TH1D *h_neg_n_StripECin2 = (TH1D*)file_in->Get("h_neg_n_StripECin2");
    h_neg_n_StripECin2->SetTitle("; Number of EC_{in} Strips");
    h_neg_n_StripECin2->Draw();
    c1->Print(Form("Figs/neg_N_ECinStrips2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/neg_N_ECinStrips2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/neg_N_ECinStrips2_%s.root", Identifier.c_str()));
    
    TH1D *h_neg_n_StripECout2 = (TH1D*)file_in->Get("h_neg_n_StripECout2");
    h_neg_n_StripECout2->SetTitle("; Number of EC_{out} Strips");
    h_neg_n_StripECout2->Draw();
    c1->Print(Form("Figs/neg_N_ECoutStrips2_%s.pdf", Identifier.c_str()));
    c1->Print(Form("Figs/neg_N_ECoutStrips2_%s.png", Identifier.c_str()));
    c1->Print(Form("Figs/neg_N_ECoutStrips2_%s.root", Identifier.c_str()));
    
    return 0;
}