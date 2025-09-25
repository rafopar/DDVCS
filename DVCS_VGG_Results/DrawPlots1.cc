/* 
 * File:   DrawPlots1.cc
 * Author: rafopar
 *
 * Created on June 9, 2022, 2:29 PM
 */

#include <cstdlib>
#include <TGraph.h>

#include <TLegend.h>
#include <TMultiGraph.h>

using namespace std;

/*
 * 
 */
void DrawPlots1() {

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.12);

    TLatex lat1;
    lat1.SetNDC();
    lat1.SetTextFont(42);
    lat1.SetTextSize(0.04);

    TGraph gr_DDVCS("results_1.dat", "%lg %*lg %*lg %*lg %*lg %lg");
    gr_DDVCS.SetMarkerColor(4);
    gr_DDVCS.SetMarkerStyle(4);
    gr_DDVCS.SetLineColor(4);

    TGraph gr_DVCS("results_DVCS_1.dat", "%lg %*lg %*lg %*lg %*lg %lg");
    gr_DVCS.SetMarkerColor(2);
    gr_DVCS.SetMarkerStyle(2);
    gr_DVCS.SetLineColor(2);

    TMultiGraph mtgr1;

    mtgr1.Add(&gr_DDVCS);
    mtgr1.Add(&gr_DVCS);
    c1->SetLogy();

    double Q2 = 2.75;
    double xB = 0.15;
    double tM = -0.3;
    double Qp2 = 3.6;
    mtgr1.Draw("AP");
    mtgr1.SetTitle("; #phi_{LH} deg; [TBD]");
    lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, t=%1.2f GeV^{2}, Q^{'2} = %1.2f GeV^{2}", Q2, xB, tM, Qp2));

    c1->Print("Figs/DDVCS_DVCS.pdf");
    c1->Print("Figs/DDVCS_DVCS.png");
    c1->Print("Figs/DDVCS_DVCS.root");

    TGraph gr_DVCS_Set2("DVCS_config2_out.dat", "%lg %*lg %*lg %*lg %*lg %lg");
    TGraph gr_DDVCS_Set2("DDVCS_config2_out.dat", "%lg %*lg %*lg %*lg %*lg %lg");
    TGraph gr_DDVCSOnly_Set2("DDVCS_Only_Setting2_out.dat", "%lg %*lg %*lg %*lg %*lg %lg");
    TGraph gr_DVCSOnly_Set2("DVCS_Only_Setting2_out.dat", "%lg %*lg %*lg %*lg %*lg %lg");

    gr_DVCS_Set2.SetMarkerColor(2);
    gr_DVCS_Set2.SetMarkerStyle(22);
    gr_DVCS_Set2.SetLineColor(2);
    gr_DVCSOnly_Set2.SetMarkerColor(2);
    gr_DVCSOnly_Set2.SetMarkerStyle(22);
    gr_DVCSOnly_Set2.SetLineColor(2);
    gr_DVCSOnly_Set2.SetLineStyle(8);

    gr_DDVCS_Set2.SetMarkerColor(4);
    gr_DDVCS_Set2.SetMarkerStyle(24);
    gr_DDVCS_Set2.SetLineColor(4);
    gr_DDVCSOnly_Set2.SetMarkerColor(4);
    gr_DDVCSOnly_Set2.SetMarkerStyle(25);
    gr_DDVCSOnly_Set2.SetLineColor(4);
    gr_DDVCSOnly_Set2.SetLineStyle(8);

    Q2 = 2.75;
    xB = 0.15;
    tM = -0.3;
    Qp2 = 0.3;


    TLegend *leg1 = new TLegend(0.4, 0.7, 0.6, 0.88);
    leg1->SetBorderSize(0);
    //leg1->AddEntry(&gr_DVCS_Set2, "DVCS");
    leg1->AddEntry(&gr_DDVCS_Set2, "DDVCS+BH");
    leg1->AddEntry(&gr_DDVCSOnly_Set2, "DDVCS");
    TMultiGraph mtgr2;
    //mtgr2.Add(&gr_DVCS_Set2);
    mtgr2.Add(&gr_DDVCS_Set2);
    mtgr2.Add(&gr_DDVCSOnly_Set2);
    mtgr2.Draw("AL");
    mtgr2.SetTitle("; #Phi_{LH} [deg]; d#sigma/dQ^{2}dt dx_{B}dQ^{'2} d#Phi [nb/GeV^{6}]");
    lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, t=%1.2f GeV^{2}, Q^{'2} = %1.2f GeV^{2}", Q2, xB, tM, Qp2));
    leg1->Draw();

    c1->Print("Figs/DDVCS_BH_Set2.pdf");
    c1->Print("Figs/DDVCS_BH_Set2.png");
    c1->Print("Figs/DDVCS_BH_Set2.root");

    TLegend *leg_DVCS = new TLegend(0.4, 0.7, 0.6, 0.88);
    leg_DVCS->SetBorderSize(0);
    leg_DVCS->AddEntry(&gr_DVCS_Set2, "DVCS + BH");
    leg_DVCS->AddEntry(&gr_DVCSOnly_Set2, "DVCS");
    TMultiGraph mtgr_DVCS;
    mtgr_DVCS.Add(&gr_DVCS_Set2);
    mtgr_DVCS.Add(&gr_DVCSOnly_Set2);
    mtgr_DVCS.Draw("AL");
    leg_DVCS->Draw();
    lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, t=%1.2f GeV^{2}", Q2, xB, tM));
    mtgr_DVCS.SetTitle("; #Phi_{LH} [deg]; d#sigma/dQ^{2}dt dx_{B} d#Phi [nb/GeV^{4}]");
    c1->Print("Figs/DVCS_BH_Set2.pdf");
    c1->Print("Figs/DVCS_BH_Set2.png");
    c1->Print("Figs/DVCS_BH_Set2.root");


    TGraph gr_DDVCSOnly_tDep_Set2("DDVCOnly_tDep_Set2_out.dat", "%*lg %lg %lg");
    TGraph gr_DDVCSPlusBH_tDep_Set2("DDVCSPlusBH_tDep_Set2_out.dat", "%*lg %lg %lg");
    TGraph gr_DVCSOnly_tDep_Set2("DVCSOnly_tDep_Set2_out.dat", "%lg %*lg %lg");
    TGraph gr_DVCSPlusBH_tDep_Set2("DVCSPlusBH_tDep_Set2_out.dat", "%lg %*lg %lg");


    gr_DVCSOnly_tDep_Set2.SetMarkerColor(2);
    gr_DVCSOnly_tDep_Set2.SetLineColor(2);
    gr_DVCSOnly_tDep_Set2.SetLineStyle(8);
    gr_DVCSPlusBH_tDep_Set2.SetMarkerColor(2);
    gr_DVCSPlusBH_tDep_Set2.SetLineColor(2);
    gr_DVCSPlusBH_tDep_Set2.SetLineStyle(1);

    TLegend *leg_DVCStDep = new TLegend(0.4, 0.7, 0.6, 0.88);
    leg_DVCStDep->SetBorderSize(0);
    leg_DVCStDep->AddEntry(&gr_DVCSPlusBH_tDep_Set2, "DVCS + BH");
    leg_DVCStDep->AddEntry(&gr_DVCSOnly_tDep_Set2, "DVCS");

    double phi_LH = 0;
    TMultiGraph *mtgr_DVCS_tDep_2 = new TMultiGraph();
    mtgr_DVCS_tDep_2->Add(&gr_DVCSOnly_tDep_Set2);
    mtgr_DVCS_tDep_2->Add(&gr_DVCSPlusBH_tDep_Set2);
    mtgr_DVCS_tDep_2->Draw("Al");
    mtgr_DVCS_tDep_2->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dt dx_{B} d#Phi [nb/GeV^{4}]");
    lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, #Phi=%1.2f deg", Q2, xB, phi_LH));
    leg_DVCStDep->Draw();
    mtgr_DVCS_tDep_2->GetXaxis()->SetLimits(0., 1.5);
    c1->Print("Figs/DVCS_tDep_Set2.pdf");
    c1->Print("Figs/DVCS_tDep_Set2.png");
    c1->Print("Figs/DVCS_tDep_Set2.root");


    gr_DDVCSOnly_tDep_Set2.SetMarkerColor(4);
    gr_DDVCSOnly_tDep_Set2.SetLineColor(4);
    gr_DDVCSOnly_tDep_Set2.SetLineStyle(8);
    gr_DDVCSPlusBH_tDep_Set2.SetMarkerColor(4);
    gr_DDVCSPlusBH_tDep_Set2.SetLineColor(4);
    gr_DDVCSPlusBH_tDep_Set2.SetLineStyle(1);

    TLegend *leg_DDVCStDep = new TLegend(0.55, 0.7, 0.75, 0.88);
    leg_DDVCStDep->SetBorderSize(0);
    leg_DDVCStDep->AddEntry(&gr_DDVCSPlusBH_tDep_Set2, "DDVCS + BH");
    leg_DDVCStDep->AddEntry(&gr_DDVCSOnly_tDep_Set2, "DDVCS");
    TMultiGraph *mtgr_DDVCS_tDep_2 = new TMultiGraph();
    mtgr_DDVCS_tDep_2->Add(&gr_DDVCSOnly_tDep_Set2);
    mtgr_DDVCS_tDep_2->Add(&gr_DDVCSPlusBH_tDep_Set2);
    mtgr_DDVCS_tDep_2->Draw("Al");
    mtgr_DDVCS_tDep_2->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dt dx_{B} d#Phi dQ^{'2} [nb/GeV^{6}]");
    lat1.DrawLatex(0.08, 0.92, Form("Q^{2}=%1.2f GeV^{2}, xB = %1.2f, Q^{'2} = %1.2f GeV^{2}, #Phi=%1.2f deg", Q2, xB, Qp2, phi_LH));
    leg_DDVCStDep->Draw();
    mtgr_DDVCS_tDep_2->GetXaxis()->SetLimits(0., 1.5);
    mtgr_DDVCS_tDep_2->SetMinimum(1.e-8);
    c1->Print("Figs/DDVCS_tDep_Set2.pdf");
    c1->Print("Figs/DDVCS_tDep_Set2.png");
    c1->Print("Figs/DDVCS_tDep_Set2.root");

    
}
