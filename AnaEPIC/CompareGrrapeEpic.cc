/* 
 * File:   CooompareGrrapeEpic.cc
 * Author: rafopar
 *
 * Created on September 23, 2023, 10:59 PM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
void CompareGrrapeEpic() {

    int runGrape = 9;
    int runEpic = 4;

    TLatex *lat1 = new TLatex();
    lat1->SetNDC();
    lat1->SetTextFont(42);
    lat1->SetTextColor(2);

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);

    TFile *file_Grape = new TFile(Form("Ana_grape_%d.root", runGrape));
    TFile *file_Epic = new TFile(Form("Ana_epic_%d.root", runEpic));

    TH1D *h_Minv2_grape = (TH1D*) file_Grape->Get("h_Minv2");
    h_Minv2_grape->SetTitle("; M(e^{-}e^{+}) [GeV]");
    TH1D *h_Minv2_epic = (TH1D*) file_Epic->Get("h_Minv2");
    h_Minv2_epic->SetLineColor(2);
    h_Minv2_epic->SetTitle("; M(e^{-}e^{+}) [GeV]");

    h_Minv2_epic->Draw("");
    h_Minv2_grape->Draw("Same");


    TRatioPlot *rp_Minv2 = new TRatioPlot(h_Minv2_epic, h_Minv2_grape);
    rp_Minv2->Draw();
    rp_Minv2->GetUpperPad()->SetBottomMargin(0);
    rp_Minv2->GetLowerPad()->SetTopMargin(0);
    rp_Minv2->GetLowerRefGraph()->SetMaximum(2.2);
    rp_Minv2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_Minv2->GetLowerRefGraph()->GetXaxis()->SetLimits(1., 3.);
    rp_Minv2->GetUpperPad()->cd();
    lat1->SetTextColor(h_Minv2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_Minv2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/Minv2Compare_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/Minv2Compare_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/Minv2Compare_grape%d_epic%d.root", runGrape, runEpic));

    c1->cd();
    TH1D *h_tM2_grape = (TH1D*) file_Grape->Get("h_tM2");
    h_tM2_grape->SetTitle("; -t [GeV^{2}]");
    TH1D *h_tM2_epic = (TH1D*) file_Epic->Get("h_tM2");
    h_tM2_epic->SetTitle("; -t [GeV^{2}]");
    h_tM2_epic->SetLineColor(2);

    TRatioPlot *rp_tM2 = new TRatioPlot(h_tM2_epic, h_tM2_grape);
    rp_tM2->Draw();
    rp_tM2->GetUpperPad()->SetBottomMargin(0);
    rp_tM2->GetLowerPad()->SetTopMargin(0);
    rp_tM2->GetLowerRefGraph()->SetMaximum(3.);
    rp_tM2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_tM2->GetLowerRefGraph()->GetXaxis()->SetLimits(0., 2.);
    rp_tM2->GetUpperPad()->cd();
    lat1->SetTextColor(h_tM2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_tM2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/tM2Compare_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/tM2Compare_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/tM2Compare_grape%d_epic%d.root", runGrape, runEpic));


    c1->cd();
    TH2D *h_Qp2_vs_Q2_2_grape = (TH2D*) file_Grape->Get("h_Qp2_vs_Q2_2");
    h_Qp2_vs_Q2_2_grape->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}]");
    h_Qp2_vs_Q2_2_grape->Draw("colz");
    c1->Print(Form("Figs/tM_grape_Run%d.pdf", runGrape));
    c1->Print(Form("Figs/tM_grape_Run%d.png", runGrape));
    c1->Print(Form("Figs/tM_grape_Run%d.root", runGrape));

    c1->cd();
    TH2D *h_Qp2_vs_Q2_2_epic = (TH2D*) file_Epic->Get("h_Qp2_vs_Q2_2");
    h_Qp2_vs_Q2_2_epic->SetTitle("; Q^{2} [GeV^{2}]; Q^{'2} [GeV^{2}]");
    h_Qp2_vs_Q2_2_epic->Draw("colz");
    c1->Print(Form("Figs/tM_epic_Run%d.pdf", runGrape));
    c1->Print(Form("Figs/tM_epic_Run%d.png", runGrape));
    c1->Print(Form("Figs/tM_epic_Run%d.root", runGrape));

    TH1D *h_Q2_2_epic = (TH1D*) h_Qp2_vs_Q2_2_epic->ProjectionX("h_Q2_2_epic", 1, h_Qp2_vs_Q2_2_epic->GetNbinsY());
    h_Q2_2_epic->SetLineColor(2);
    TH1D *h_Q2_2_grape = (TH1D*) h_Qp2_vs_Q2_2_grape->ProjectionX("h_Q2_2_grape", 1, h_Qp2_vs_Q2_2_grape->GetNbinsY());

    TRatioPlot *rp_Q2_2 = new TRatioPlot(h_Q2_2_epic, h_Q2_2_grape);
    rp_Q2_2->Draw();
    rp_Q2_2->GetUpperPad()->SetBottomMargin(0);
    rp_Q2_2->GetLowerPad()->SetTopMargin(0);
    rp_Q2_2->GetLowerRefGraph()->SetMaximum(3.);
    rp_Q2_2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_Q2_2->GetLowerRefGraph()->GetXaxis()->SetLimits(0., 5.);
    rp_Q2_2->GetUpperPad()->cd();
    lat1->SetTextColor(h_Q2_2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_Q2_2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/rp_Q2_2_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/rp_Q2_2_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/rp_Q2_2_grape%d_epic%d.root", runGrape, runEpic));

    c1->cd();
    TH2D *h_th_P_em2_grape = (TH2D*) file_Grape->Get("h_th_P_em2");
    h_th_P_em2_grape->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-)} [deg]");
    h_th_P_em2_grape->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_grape_Run_%d.pdf", runGrape));
    c1->Print(Form("Figs/th_P_em2_grape_Run_%d.png", runGrape));
    c1->Print(Form("Figs/th_P_em2_grape_Run_%d.root", runGrape));

    TH2D *h_th_P_em2_epic = (TH2D*) file_Epic->Get("h_th_P_em2");
    h_th_P_em2_epic->SetTitle("; P(e^{-}) [GeV]; #theta(e^{-)} [deg]");
    h_th_P_em2_epic->Draw("colz");
    c1->Print(Form("Figs/th_P_em2_epic_Run_%d.pdf", runEpic));
    c1->Print(Form("Figs/th_P_em2_epic_Run_%d.png", runEpic));
    c1->Print(Form("Figs/th_P_em2_epic_Run_%d.root", runEpic));

    TH2D *h_Q2_xB2_grape = (TH2D*) file_Grape->Get("h_Q2_xB2");
    h_Q2_xB2_grape->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2_grape->SetMarkerColor(4);

    TH2D *h_Q2_xB2_epic = (TH2D*) file_Epic->Get("h_Q2_xB2");
    h_Q2_xB2_epic->SetTitle("; x_{B}; Q^{2} [GeV^{2}]");
    h_Q2_xB2_epic->SetMarkerColor(2);

    h_Q2_xB2_grape->Draw("colz");
    c1->Print(Form("Figs/Q2_xB2_grape_Run_%d.pdf", runGrape));
    c1->Print(Form("Figs/Q2_xB2_grape_Run_%d.png", runGrape));
    c1->Print(Form("Figs/Q2_xB2_grape_Run_%d.root", runGrape));

    h_Q2_xB2_epic->Draw("colz");
    c1->Print(Form("Figs/Q2_xB2_epic_Run_%d.pdf", runEpic));
    c1->Print(Form("Figs/Q2_xB2_epic_Run_%d.png", runEpic));
    c1->Print(Form("Figs/Q2_xB2_epic_Run_%d.root", runEpic));

    TH1D *h_xB2_grape = h_Q2_xB2_grape->ProjectionX("h_xB2_grape", 1, h_Q2_xB2_grape->GetNbinsY());
    h_xB2_grape->SetLineColor(4);
    TH1D *h_xB2_epic = h_Q2_xB2_epic->ProjectionX("h_xB2_epic", 1, h_Q2_xB2_epic->GetNbinsY());
    h_xB2_grape->SetLineColor(2);

    TRatioPlot *rp_xB2 = new TRatioPlot(h_xB2_epic, h_xB2_grape);
    rp_xB2->Draw();
    rp_xB2->Draw();
    rp_xB2->GetUpperPad()->SetBottomMargin(0);
    rp_xB2->GetLowerPad()->SetTopMargin(0);
    rp_xB2->GetLowerRefGraph()->SetMaximum(3.);
    rp_xB2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_xB2->GetLowerRefGraph()->GetXaxis()->SetLimits(0., 0.4);
    rp_xB2->GetUpperPad()->cd();
    lat1->SetTextColor(h_xB2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_xB2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/xB2_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/xB2_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/xB2_grape%d_epic%d.root", runGrape, runEpic));

    c1->cd();
    TH2D *h_th_P_mum2_grape = (TH2D*) file_Grape->Get("h_th_P_mum2");
    h_th_P_mum2_grape->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-)} [deg]");
    h_th_P_mum2_grape->Draw("colz");
    c1->Print(Form("Figs/th_P_mum2_grape_Run_%d.pdf", runGrape));
    c1->Print(Form("Figs/th_P_mum2_grape_Run_%d.png", runGrape));
    c1->Print(Form("Figs/th_P_mum2_grape_Run_%d.root", runGrape));

    TH2D *h_th_P_mum2_epic = (TH2D*) file_Epic->Get("h_th_P_mum2");
    h_th_P_mum2_epic->SetTitle("; P(#mu^{-}) [GeV]; #theta(#mu^{-)} [deg]");
    h_th_P_mum2_epic->Draw("colz");
    c1->Print(Form("Figs/th_P_mum2_epic_Run_%d.pdf", runEpic));
    c1->Print(Form("Figs/th_P_mum2_epic_Run_%d.png", runEpic));
    c1->Print(Form("Figs/th_P_mum2_epic_Run_%d.root", runEpic));

    TH1D *h_P_mum2_grape = (TH1D*) h_th_P_mum2_grape->ProjectionX("h_P_mum2_grape", 1, h_th_P_mum2_grape->GetNbinsY());
    h_P_mum2_grape->SetLineColor(4);
    TH1D *h_P_mum2_epic = (TH1D*) h_th_P_mum2_epic->ProjectionX("h_P_mum2_epic", 1, h_th_P_mum2_epic->GetNbinsY());
    h_P_mum2_epic->SetLineColor(2);

    TRatioPlot *rp_P_mum2 = new TRatioPlot(h_P_mum2_epic, h_P_mum2_grape);
    rp_P_mum2->Draw();
    rp_P_mum2->GetUpperPad()->SetBottomMargin(0);
    rp_P_mum2->GetLowerPad()->SetTopMargin(0);
    rp_P_mum2->GetLowerRefGraph()->SetMaximum(4.);
    rp_P_mum2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_P_mum2->GetLowerRefGraph()->GetXaxis()->SetLimits(0.9, 8.);
    rp_P_mum2->GetUpperPad()->cd();
    lat1->SetTextColor(h_P_mum2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_P_mum2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/P_mum2_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/P_mum2_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/P_mum2_grape%d_epic%d.root", runGrape, runEpic));

    c1->cd();
    TH2D *h_th_P_mup2_grape = (TH2D*) file_Grape->Get("h_th_P_mup2");
    h_th_P_mup2_grape->SetTitle("; P(#mu^{+}) [GeV]; #theta(#mu^{+)} [deg]");
    h_th_P_mup2_grape->Draw("colz");
    c1->Print(Form("Figs/th_P_mup2_grape_Run_%d.pdf", runGrape));
    c1->Print(Form("Figs/th_P_mup2_grape_Run_%d.png", runGrape));
    c1->Print(Form("Figs/th_P_mup2_grape_Run_%d.root", runGrape));

    TH2D *h_th_P_mup2_epic = (TH2D*) file_Epic->Get("h_th_P_mup2");
    h_th_P_mup2_epic->SetTitle("; P(#mu^{+}) [GeV]; #theta(#mu^{+)} [deg]");
    h_th_P_mup2_epic->Draw("colz");
    c1->Print(Form("Figs/th_P_mup2_epic_Run_%d.pdf", runEpic));
    c1->Print(Form("Figs/th_P_mup2_epic_Run_%d.png", runEpic));
    c1->Print(Form("Figs/th_P_mup2_epic_Run_%d.root", runEpic));

    TH1D *h_P_mup2_grape = (TH1D*) h_th_P_mup2_grape->ProjectionX("h_P_mup2_grape", 1, h_th_P_mup2_grape->GetNbinsY());
    h_P_mup2_grape->SetLineColor(4);
    TH1D *h_P_mup2_epic = (TH1D*) h_th_P_mup2_epic->ProjectionX("h_P_mup2_epic", 1, h_th_P_mup2_epic->GetNbinsY());
    h_P_mup2_epic->SetLineColor(2);

    TRatioPlot *rp_P_mup2 = new TRatioPlot(h_P_mup2_epic, h_P_mup2_grape);
    
    rp_P_mup2->Draw();
    rp_P_mup2->GetUpperPad()->SetBottomMargin(0);
    rp_P_mup2->GetLowerPad()->SetTopMargin(0);
    rp_P_mup2->GetLowerRefGraph()->SetMaximum(4.);
    rp_P_mup2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_P_mup2->GetLowerRefGraph()->GetXaxis()->SetLimits(0.9, 8.);
    rp_P_mup2->GetUpperPad()->cd();
    lat1->SetTextColor(h_P_mup2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_P_mup2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/P_mup2_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/P_mup2_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/P_mup2_grape%d_epic%d.root", runGrape, runEpic));

    c1->Clear();
    TH2D *h_th_P_p2_grape = (TH2D*) file_Grape->Get("h_th_P_p2");
    h_th_P_p2_grape->SetTitle("; P(proton) [GeV]; #theta(#mu^{-)} [deg]");
    h_th_P_p2_grape->Draw("colz");
    c1->Print(Form("Figs/th_P_p2_grape_Run_%d.pdf", runGrape));
    c1->Print(Form("Figs/th_P_p2_grape_Run_%d.png", runGrape));
    c1->Print(Form("Figs/th_P_p2_grape_Run_%d.root", runGrape));

    
    TH2D *h_th_P_p2_epic = (TH2D*) file_Epic->Get("h_th_P_p2");
    h_th_P_p2_epic->SetTitle("; P(proton) [GeV]; #theta(#mu^{-)} [deg]");
    h_th_P_p2_epic->Draw("colz");
    c1->Print(Form("Figs/th_P_p2_epic_Run_%d.pdf", runEpic));
    c1->Print(Form("Figs/th_P_p2_epic_Run_%d.png", runEpic));
    c1->Print(Form("Figs/th_P_p2_epic_Run_%d.root", runEpic));

    TH1D *h_P_p2_grape = (TH1D*) h_th_P_p2_grape->ProjectionX("h_P_p2_grape", 1, h_th_P_p2_grape->GetNbinsY());
    h_P_p2_grape->SetLineColor(4);
    TH1D *h_P_p2_epic = (TH1D*) h_th_P_p2_epic->ProjectionX("h_P_p2_epic", 1, h_th_P_p2_epic->GetNbinsY());
    h_P_p2_epic->SetLineColor(2);

    TRatioPlot *rp_P_p2 = new TRatioPlot(h_P_p2_epic, h_P_p2_grape);
    
    rp_P_p2->Draw();
    rp_P_p2->GetUpperPad()->SetBottomMargin(0);
    rp_P_p2->GetLowerPad()->SetTopMargin(0);
    rp_P_p2->GetLowerRefGraph()->SetMaximum(3.);
    rp_P_p2->GetLowerRefGraph()->SetMinimum(0.9);
    rp_P_p2->GetLowerRefGraph()->GetXaxis()->SetLimits(0.1, 2.);
    rp_P_p2->GetUpperPad()->cd();
    lat1->SetTextColor(h_P_p2_epic->GetLineColor());
    lat1->DrawLatex(0.55, 0.8, "EPIC");
    lat1->SetTextColor(h_P_p2_grape->GetLineColor());
    lat1->DrawLatex(0.55, 0.7, "GRAPE");
    c1->Print(Form("Figs/P_p2_grape%d_epic%d.pdf", runGrape, runEpic));
    c1->Print(Form("Figs/P_p2_grape%d_epic%d.png", runGrape, runEpic));
    c1->Print(Form("Figs/P_p2_grape%d_epic%d.root", runGrape, runEpic));

}
