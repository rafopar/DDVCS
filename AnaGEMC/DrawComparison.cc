/* 
 * File:   DrawComparison.cc
 * Author: rafopar
 *
 * Created on November 25, 2024, 10:51 AM
 */

#include <map>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLine.h>
#include <cstdlib>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

template <typename T> void FormatHist(T h);
/*
 * 
 */
int main(int argc, char** argv) {

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.05);
    c1->SetRightMargin(0.02);
    TCanvas *c2 = new TCanvas("c2", "", 950, 950);
    c2->SetTopMargin(0.02);
    c2->SetRightMargin(0.02);
    TCanvas *c3 = new TCanvas("c3", "", 950, 950);
    c3->SetTopMargin(0.02);
    c3->SetRightMargin(0.02);
    TCanvas *c4 = new TCanvas("c4", "", 950, 950);
    c4->SetTopMargin(0.02);
    c4->SetRightMargin(0.02);

    //std::vector<int> v_runs = {17, 19, 9113, 9114};
    std::vector<int> v_runs = {17, 19};
    std::vector<int> v_22GeVruns = {18, 20, 9115, 9116};

    std::map<int, TH1D*> m_h_MX2;
    std::map<int, TH1D*> m_h_MX2_2;
    std::map<int, TH1D*> m_h_MX2_3;
    std::map<int, TH1D*> m_h_Mmumu2;

    std::map<int, TH1D*> m_h_MX2_22GeV;
    std::map<int, TH1D*> m_h_MX2_2_22GeV;
    std::map<int, TH1D*> m_h_MX2_3_22GeV;
    std::map<int, TH1D*> m_h_Mmumu2_22GeV;

    std::map<int, int> m_colHist;
    std::map<int, int> m_colHist_22GeV;

    m_colHist[17] = 4;
    m_colHist[19] = 2;
    m_colHist[9113] = 6;
    m_colHist[9114] = 8;

    m_colHist_22GeV[18] = 4;
    m_colHist_22GeV[20] = 2;
    m_colHist_22GeV[9115] = 6;
    m_colHist_22GeV[9116] = 8;

    std::map<int, std::string> m_Description;
    m_Description[17] = "Elastic";
    m_Description[19] = "Quasi-elastic DDVCS";
    m_Description[9113] = "#rho #rightarrow #mu^{-}#mu^{+}";
    m_Description[9114] = "#rho #rightarrow #pi^{-}#pi^{+}";

    std::map<int, std::string> m_Description_22GeV;
    m_Description_22GeV[18] = "Elastic";
    m_Description_22GeV[20] = "Quasi-elastic DDVCS";
    m_Description_22GeV[9115] = "#rho #rightarrow #mu^{-}#mu^{+}";
    m_Description_22GeV[9116] = "#rho #rightarrow #pi^{-}#pi^{+}";

    const double Mx2_Max = 1.5;
    const double Mx2_Min = 0.4;

    TLine *line1 = new TLine();

    double max_MX2 = 0;
    TLegend *leg_MX2 = new TLegend(0.5, 0.6, 0.97, 0.94);
    leg_MX2->SetBorderSize(0);
    TLegend *leg_Mmumu2 = new TLegend(0.57, 0.6, 0.97, 0.94);
    leg_Mmumu2->SetBorderSize(0);
    for (int i = 0; i < v_runs.size(); i++) {

        int run = v_runs.at(i);

        TFile *file_in = new TFile(Form("AnaDDVCS_Run_%d.root", run));

        c1->cd();
        m_h_MX2[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr1");
        m_h_MX2[run]->SetTitle("; M_{X}^{2} [GeV^{2}]");
        m_h_MX2[run]->SetLineWidth(2);
        m_h_MX2[run]->SetLineColor(m_colHist[run]);

        leg_MX2->AddEntry(m_h_MX2[run], m_Description[run].c_str());
        m_h_MX2[run]->Draw("hist Same");

        c2->cd();
        m_h_Mmumu2[run] = (TH1D*) file_in->Get("h_Mumu_Corr_AngleFix2");
        m_h_Mmumu2[run]->SetTitle("; M(#mu^{-}#mu^{+}) [GeV]");
        m_h_Mmumu2[run]->SetLineWidth(2);
        m_h_Mmumu2[run]->SetLineColor(m_colHist[run]);
        leg_Mmumu2->AddEntry(m_h_Mmumu2[run], m_Description[run].c_str());
        m_h_Mmumu2[run]->Draw("hist same");

        c3->cd();
        m_h_MX2_2[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr2");
        m_h_MX2_2[run]->SetTitle("; M_{X}^{2} [GeV^{2}]");
        m_h_MX2_2[run]->SetLineWidth(2);
        m_h_MX2_2[run]->SetLineColor(m_colHist[run]);
        m_h_MX2_2[run]->Draw("hist Same");

        c4->cd();
        m_h_MX2_3[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr3");
        m_h_MX2_3[run]->SetTitle("; M_{X}^{2} [GeV^{2}]");
        m_h_MX2_3[run]->SetLineWidth(2);
        m_h_MX2_3[run]->SetLineColor(m_colHist[run]);
        m_h_MX2_3[run]->Draw("hist Same");
    }


    //m_h_MX2[ v_runs.at(0) ]->SetMaximum(1.05*max_MX2);
    c1->cd();
    line1->SetLineWidth(2);
    line1->SetLineStyle(9);
    line1->DrawLine(Mx2_Min, 0, Mx2_Min, m_h_MX2[ v_runs.at(0) ]->GetMaximum());
    line1->DrawLine(Mx2_Max, 0, Mx2_Max, m_h_MX2[ v_runs.at(0) ]->GetMaximum());
    leg_MX2->Draw();
    c1->Print("Figs/DDVCS_Mx2_comparisons.pdf");
    c1->Print("Figs/DDVCS_Mx2_comparisons.png");
    c1->Print("Figs/DDVCS_Mx2_comparisons.root");

    c3->cd();
    line1->DrawLine(Mx2_Min, 0, Mx2_Min, m_h_MX2_2[ v_runs.at(0) ]->GetMaximum());
    line1->DrawLine(Mx2_Max, 0, Mx2_Max, m_h_MX2_2[ v_runs.at(0) ]->GetMaximum());
    leg_MX2->Draw();
    c3->Print("Figs/DDVCS_Mx2_2_comparisons.pdf");
    c3->Print("Figs/DDVCS_Mx2_2_comparisons.png");
    c3->Print("Figs/DDVCS_Mx2_2_comparisons.root");

    c4->cd();
    line1->DrawLine(Mx2_Min, 0, Mx2_Min, m_h_MX2_3[ v_runs.at(0) ]->GetMaximum());
    line1->DrawLine(Mx2_Max, 0, Mx2_Max, m_h_MX2_3[ v_runs.at(0) ]->GetMaximum());
    leg_MX2->Draw();
    c4->Print("Figs/DDVCS_Mx2_3_comparisons.pdf");
    c4->Print("Figs/DDVCS_Mx2_3_comparisons.png");
    c4->Print("Figs/DDVCS_Mx2_3_comparisons.root");

    c2->cd();
    leg_Mmumu2->Draw();
    c2->Print("Figs/DDVCS_Mmumu_comparisons.pdf");
    c2->Print("Figs/DDVCS_Mmumu_comparisons.png");
    c2->Print("Figs/DDVCS_Mmumu_comparisons.root");



    c1->Clear();
    c2->Clear();
    c3->Clear();
    c4->Clear();
    leg_MX2->Clear();
    leg_Mmumu2->Clear();

    for (int i = 0; i < v_22GeVruns.size(); i++) {

        int run = v_22GeVruns.at(i);
        TFile *file_in = new TFile(Form("AnaDDVCS_Run_%d.root", run));

        c1->cd();
        m_h_MX2_22GeV[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr1");
        m_h_MX2_22GeV[run]->SetTitle("; M_{X}^{2} [GeV^{2}]");
        m_h_MX2_22GeV[run]->SetLineWidth(2);
        m_h_MX2_22GeV[run]->SetLineColor(m_colHist_22GeV[run]);

        leg_MX2->AddEntry(m_h_MX2_22GeV[run], m_Description_22GeV[run].c_str());
        if (i == 0) {
            m_h_MX2_22GeV[run]->Draw("hist");
        } else {
            m_h_MX2_22GeV[run]->Draw("hist Same");
        }

        c2->cd();
        m_h_Mmumu2_22GeV[run] = (TH1D*) file_in->Get("h_Mumu_Corr_AngleFix2");
        m_h_Mmumu2_22GeV[run]->SetTitle("; M(#mu^{-}#mu^{+}) [GeV]");
        m_h_Mmumu2_22GeV[run]->SetLineWidth(2);
        m_h_Mmumu2_22GeV[run]->SetLineColor(m_colHist_22GeV[run]);
        leg_Mmumu2->AddEntry(m_h_Mmumu2_22GeV[run], m_Description_22GeV[run].c_str());
        if (i == 0) {
            m_h_Mmumu2_22GeV[run]->Draw("hist");
        } else {
            m_h_Mmumu2_22GeV[run]->Draw("hist same");
        }

        c3->cd();
        m_h_MX2_2_22GeV[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr2");
        m_h_MX2_2_22GeV[run]->SetTitle("; M_{X}^{2} [GeV^{2}]");
        m_h_MX2_2_22GeV[run]->SetLineWidth(2);
        m_h_MX2_2_22GeV[run]->SetLineColor(m_colHist_22GeV[run]);
        m_h_MX2_2_22GeV[run]->Draw("hist Same");

        c4->cd();
        m_h_MX2_3_22GeV[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr3");
        m_h_MX2_3_22GeV[run]->SetTitle("; M_{X}^{2} [GeV^{2}]");
        m_h_MX2_3_22GeV[run]->SetLineWidth(2);
        m_h_MX2_3_22GeV[run]->SetLineColor(m_colHist_22GeV[run]);
        m_h_MX2_3_22GeV[run]->Draw("hist Same");
    }


    //m_h_MX2[ v_runs.at(0) ]->SetMaximum(1.05*max_MX2);
    c1->cd();
    line1->SetLineWidth(2);
    line1->SetLineStyle(9);
    line1->DrawLine(Mx2_Min, 0, Mx2_Min, m_h_MX2_22GeV[ v_22GeVruns.at(0) ]->GetMaximum());
    line1->DrawLine(Mx2_Max, 0, Mx2_Max, m_h_MX2_22GeV[ v_22GeVruns.at(0) ]->GetMaximum());

    leg_MX2->Draw();
    c1->Print("Figs/DDVCS_Mx2_comparisons_22GeV.pdf");
    c1->Print("Figs/DDVCS_Mx2_comparisons_22GeV.png");
    c1->Print("Figs/DDVCS_Mx2_comparisons_22GeV.root");

    c3->cd();
    line1->DrawLine(Mx2_Min, 0, Mx2_Min, m_h_MX2_2_22GeV[ v_22GeVruns.at(0) ]->GetMaximum());
    line1->DrawLine(Mx2_Max, 0, Mx2_Max, m_h_MX2_2_22GeV[ v_22GeVruns.at(0) ]->GetMaximum());
    leg_MX2->Draw();
    c3->Print("Figs/DDVCS_Mx2_2_comparisons_22GeV.pdf");
    c3->Print("Figs/DDVCS_Mx2_2_comparisons_22GeV.png");
    c3->Print("Figs/DDVCS_Mx2_2_comparisons_22GeV.root");

    c4->cd();
    line1->DrawLine(Mx2_Min, 0, Mx2_Min, m_h_MX2_3_22GeV[ v_22GeVruns.at(0) ]->GetMaximum());
    line1->DrawLine(Mx2_Max, 0, Mx2_Max, m_h_MX2_3_22GeV[ v_22GeVruns.at(0) ]->GetMaximum());
    leg_MX2->Draw();
    c4->Print("Figs/DDVCS_Mx2_3_comparisons_22GeV.pdf");
    c4->Print("Figs/DDVCS_Mx2_3_comparisons_22GeV.png");
    c4->Print("Figs/DDVCS_Mx2_3_comparisons_22GeV.root");

    c2->cd();
    leg_Mmumu2->Draw();
    c2->Print("Figs/DDVCS_Mmumu_comparisons_22GeV.pdf");
    c2->Print("Figs/DDVCS_Mmumu_comparisons_22GeV.png");
    c2->Print("Figs/DDVCS_Mmumu_comparisons_22GeV.root");

    
    /*
     * Plot Q'2 distributions from 22 GeV and 10.6 GeV on top of each other
     */
    
    TFile *file_10p6GeV = new TFile("AnaDDVCS_Run_17.root");
    TFile *file_22GeV = new TFile("AnaDDVCS_Run_18.root");
    
    int n_x_xi_bins = 2;
    
    c1->cd();
    c1->Clear();
    for( int i = 0; i < n_x_xi_bins; i++ ){
        TH1D *h_M_mumu_10p6GeV = (TH1D*)file_10p6GeV->Get(Form("h_M_mumu_%d", i));
        h_M_mumu_10p6GeV->SetTitle("; M(#mu^{-}#mu^{+}) [GeV]");
        h_M_mumu_10p6GeV->SetLineColor(2);
        h_M_mumu_10p6GeV->SetLineWidth(3);
        FormatHist(h_M_mumu_10p6GeV);
        
        TH1D *h_M_mumu_22GeV = (TH1D*)file_22GeV->Get(Form("h_M_mumu_%d", i));
        h_M_mumu_22GeV->SetTitle("; M(#mu^{-}#mu^{+}) [GeV]");
        h_M_mumu_22GeV->SetLineColor(4);
        h_M_mumu_22GeV->SetLineWidth(4);
        FormatHist(h_M_mumu_22GeV);
        
        h_M_mumu_10p6GeV->SetMaximum(1.05*TMath::Max(h_M_mumu_10p6GeV->GetMaximum(), h_M_mumu_22GeV->GetMaximum()) );
        h_M_mumu_10p6GeV->Draw("hist");
        h_M_mumu_22GeV->Draw("hist same");
        c1->Print(Form("Figs/M_mumu_Comparisons_10p6AND22GeV_Bin_%d.pdf", i));
        c1->Print(Form("Figs/M_mumu_Comparisons_10p6AND22GeV_Bin_%d.png", i));
        c1->Print(Form("Figs/M_mumu_Comparisons_10p6AND22GeV_Bin_%d.root", i));
    }
    
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