/* 
 * File:   DrawAngularResolutions.cc
 * Author: rafopar
 *
 * Created on July 3, 2025, 5:26 PM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
void DrawAngularResolutions() {

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    
    TGraph *gr_sigm_th_mum;
    TGraph *gr_sigm_th_mup;
    TGraph *gr_sigm_phi_mum;
    TGraph *gr_sigm_phi_mup;
    
    gr_sigm_th_mum = new TGraph("AngularResolutions_60cm_Lead.dat", "%lg %lg");
    gr_sigm_th_mum->SetMarkerStyle(20);
    gr_sigm_th_mum->SetMarkerSize(2);
    gr_sigm_th_mum->SetMarkerColor(2);
    gr_sigm_th_mum->SetTitle("; P [GeV]; #sigma (#theta) [deg]");
    gr_sigm_th_mum->Draw("AP");

    c1->Clear();
    gr_sigm_th_mup = new TGraph("AngularResolutions_60cm_Lead.dat", "%lg %*lg %lg");
    gr_sigm_th_mup->SetMarkerStyle(21);
    gr_sigm_th_mup->SetMarkerSize(2);
    gr_sigm_th_mup->SetMarkerColor(4);
    gr_sigm_th_mup->SetTitle("; P [GeV]; #sigma (#theta) [deg]");
    gr_sigm_th_mup->Draw("AP");
    
    TLegend *leg_th = new TLegend(0.55, 0.65, 0.9, 0.9);
    leg_th->SetBorderSize(0);
    leg_th->AddEntry(gr_sigm_th_mum, "#mu^{-}");
    leg_th->AddEntry(gr_sigm_th_mup, "#mu^{+}");
    
    c1->Clear();
    TMultiGraph *mtgr_th = new TMultiGraph();
    mtgr_th->Add(gr_sigm_th_mum);
    mtgr_th->Add(gr_sigm_th_mup);
    mtgr_th->Draw("AP");
    //mtgr_th->SetMinimum(0);
    mtgr_th->SetTitle("; P [GeV]; #sigma (#theta) [deg]");
    leg_th->Draw();
    c1->Print("Figs/th_resolutions.pdf");
    c1->Print("Figs/th_resolutions.png");
    c1->Print("Figs/th_resolutions.root");
    
    gr_sigm_phi_mum = new TGraph("AngularResolutions_60cm_Lead.dat", "%lg %*lg %*lg %lg");
    gr_sigm_phi_mum->SetMarkerStyle(20);
    gr_sigm_phi_mum->SetMarkerSize(2);
    gr_sigm_phi_mum->SetMarkerColor(2);
    gr_sigm_phi_mum->SetTitle("; P [GeV]; #sigma (#phi) [deg]");
    gr_sigm_phi_mum->Draw("AP");

    c1->Clear();
    gr_sigm_phi_mup = new TGraph("AngularResolutions_60cm_Lead.dat", "%lg %*lg %*lg %*lg %lg");
    gr_sigm_phi_mup->SetMarkerStyle(21);
    gr_sigm_phi_mup->SetMarkerSize(2);
    gr_sigm_phi_mup->SetMarkerColor(4);
    gr_sigm_phi_mup->SetTitle("; P [GeV]; #sigma (#phi) [deg]");
    gr_sigm_phi_mup->Draw("AP");
    
    TLegend *leg_phi = new TLegend(0.55, 0.65, 0.9, 0.9);
    leg_phi->SetBorderSize(0);
    leg_phi->AddEntry(gr_sigm_phi_mum, "#mu^{-}");
    leg_phi->AddEntry(gr_sigm_phi_mup, "#mu^{+}");
    
    c1->Clear();
    TMultiGraph *mtgr_phi = new TMultiGraph();
    mtgr_phi->Add(gr_sigm_phi_mum);
    mtgr_phi->Add(gr_sigm_phi_mup);
    mtgr_phi->Draw("AP");
    //mtgr_phi->SetMinimum(0);
    mtgr_phi->SetTitle("; P [GeV]; #sigma (#phi) [deg]");
    leg_phi->Draw();
    c1->Print("Figs/phi_resolutions.pdf");
    c1->Print("Figs/phi_resolutions.png");
    c1->Print("Figs/phi_resolutions.root");
}

