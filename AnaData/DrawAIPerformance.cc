//
// Created by rafopar on 1/23/26.
//

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>

#include <iostream>

using namespace std;

TH1D* GetContaminationHistogram(TH1D* h_True, TH1D* h_False);

int main (int argc, char *argv[]) {

    if (argc != 2) {
        cout << "Usage: DrawAIPerformance.exe <RUN>" << argv[0] << endl;
    }

    int run = atoi(argv[1]);

    auto c1 = new TCanvas("c1", "", 950, 950);
    c1->SetRightMargin(0.02);

    auto file_in = new TFile(Form("HistData/Hists_DDVCS_Run_%d.root", run));


    auto h_Score_True1 = dynamic_cast<TH1D*>(file_in->Get("h_Score_True1"));
    h_Score_True1->SetTitle(";Score");
    h_Score_True1->SetLineColor(kBlue);
    h_Score_True1->SetLineWidth(2);
    h_Score_True1->Draw("hist");

    auto h_Score_False1 = dynamic_cast<TH1D*>(file_in->Get("h_Score_False1"));
    h_Score_False1->SetTitle(";Score");
    h_Score_False1->SetLineColor(kRed);
    h_Score_False1->SetLineWidth(2);
    h_Score_False1->Draw("hist Same");
    c1->Print(Form("Figs/score_TrueFalse_Run_%d.pdf", run));
    c1->Print(Form("Figs/score_TrueFalse_Run_%d.png", run));
    c1->Print(Form("Figs/score_TrueFalse_Run_%d.root", run));

    auto h_Contamination = GetContaminationHistogram(h_Score_True1, h_Score_False1);
    h_Contamination->Draw("hist");
    c1->Print(Form("Figs/Containation_Vs_Score_Run_%d.pdf", run));
    c1->Print(Form("Figs/Containation_Vs_Score_Run_%d.png", run));
    c1->Print(Form("Figs/Containation_Vs_Score_Run_%d.root", run));

    return 0;
}

TH1D* GetContaminationHistogram(TH1D* h_True, TH1D* h_False) {
    double frac = h_False->Integral()/h_True->Integral();
    double NTot = h_False->Integral() + h_True->Integral();

    auto h_Contamination = dynamic_cast<TH1D*>(h_True->Clone("h_Contamination"));
    h_Contamination->Reset();

    for (auto i = 0; i < h_True->GetNbinsX(); i++) {

        double N_WrongChoice_True = h_True->Integral(i+1, h_True->GetNbinsX());
        double N_RightChoice_True = h_True->Integral(1, i+1);

        double N_WrongChoice_False = h_False->Integral(1, i+1);
        double N_RightChoice_False = h_False->Integral(i+1, h_False->GetNbinsX());

        double Contamination = (N_WrongChoice_False + N_WrongChoice_True)/NTot;
        h_Contamination->SetBinContent(i+1, Contamination);
    }
    return h_Contamination;
}