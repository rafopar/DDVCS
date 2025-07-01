/* 
 * File:   DrawInclusiveCoincidence.cc
 * Author: rafopar
 *
 * Created on March 3, 2025, 8:28 PM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
void DrawInclusiveCoincidence() {

    gStyle->SetOptStat(0);

    const int RunBH = 17;
    const int RunElasticMerge = 2023;
    const int RunQuasiElasticMerge = 1022;

    const double MM2_max = 1.5; // From AnaGrape
    const double MM2_min = 0.4; // From AnaGrape


    const int nRuns = 3;

    int runs[nRuns] = {RunBH, RunElasticMerge, RunQuasiElasticMerge};

    std::map<int, TH1D*> mh_Mmis_;

    for (int i = 0; i < nRuns; i++) {

        int run = runs[i];
        TFile *file_in = new TFile(Form("AnaDDVCS_Run_%d.root", run), "Read");
        mh_Mmis_[run] = (TH1D*) file_in->Get("h_Mmis_AngleFixcorr1");
        mh_Mmis_[run]->SetLineWidth(3);
    }

    TH1D *h_BgrSum = (TH1D*) mh_Mmis_[RunElasticMerge]->Clone("h_BgrSum");
    h_BgrSum->Add(mh_Mmis_[RunQuasiElasticMerge]);
    h_BgrSum->SetLineColor(96);

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);

    mh_Mmis_[RunBH]->SetTitle("; M_{X}^{2} [GeV^{2}]");
    mh_Mmis_[RunBH]->SetLineColor(4);
    mh_Mmis_[RunElasticMerge]->SetLineColor(2);
    mh_Mmis_[RunQuasiElasticMerge]->SetLineColor(6);

    TLegend *leg1 = new TLegend(0.4, 0.65, 0.95, 0.8);
    leg1->SetBorderSize(0);
    leg1->AddEntry(mh_Mmis_[RunBH], "BH");
    //leg1->AddEntry(mh_Mmis_[RunElasticMerge], "Coincidence w/ elast");
    //leg1->AddEntry(mh_Mmis_[RunQuasiElasticMerge], "Coincidence w/ quasi-elast");
    leg1->AddEntry(h_BgrSum, "Coincidence w/ el + quasi-elast");

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextColor(96);
    lat1->SetTextFont(42);
    
    int binLeft1 = mh_Mmis_[RunBH]->FindBin(MM2_min);
    int binRight1 = mh_Mmis_[RunBH]->FindBin(MM2_max);
    
    double integ_BH = mh_Mmis_[RunBH]->Integral(binLeft1, binRight1);
    double integ_BgrSum = h_BgrSum->Integral(binLeft1, binRight1);
    
    mh_Mmis_[RunBH]->Draw();
    //mh_Mmis_[RunElasticMerge]->Draw("Same");
    //mh_Mmis_[RunQuasiElasticMerge]->Draw("Same");
    h_BgrSum->Draw("Same");
    leg1->Draw();
    lat1->DrawLatex(0.45, 0.93, Form("Tot contribution is %1.1f %%", 100.*integ_BgrSum/integ_BH ));
    c1->Print(Form("Figs/Inclusive_Coincidence_LinScale.pdf"));
    c1->Print(Form("Figs/Inclusive_Coincidence_LinScale.png"));
    c1->Print(Form("Figs/Inclusive_Coincidence_LinScale.root"));

    c1->SetLogy();
    c1->Print(Form("Figs/Inclusive_Coincidence_LogScale.pdf"));
    c1->Print(Form("Figs/Inclusive_Coincidence_LogScale.png"));
    c1->Print(Form("Figs/Inclusive_Coincidence_LogScale.root"));


}