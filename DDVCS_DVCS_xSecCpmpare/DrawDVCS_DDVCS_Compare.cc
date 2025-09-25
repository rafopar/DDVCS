/* 
 * File:   DrawDVCS_DDVCS_Compare.cc
 * Author: rafopar
 *
 * Created on February 6, 2025, 6:20 PM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
void DrawDVCS_DDVCS_Compare() {

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.14);
    c1->SetLogy();
    
    const double Eb = 10.6;
    const double Q2 = 2.75;
    const double xB = 0.15;
    const double Qp2 = 1.4;


    TGraph *gr_DDVCS_tDep = new TGraph("DDVCS_10p6_xsec.dat", "%*lg %lg %lg");

    gr_DDVCS_tDep->SetLineColor(4);
    gr_DDVCS_tDep->SetLineStyle(8);
    gr_DDVCS_tDep->SetLineWidth(4);
    gr_DDVCS_tDep->Draw("AL");
    gr_DDVCS_tDep->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dtdx_{B}d#Phi dQ^{'2}[nb/GeV^{6}]");

    TGraph *gr_DDVCSPlusBH_tDep = new TGraph("DDVCSPlusBH_10p6GeV_xsec.dat", "%*lg %lg %lg");

    gr_DDVCSPlusBH_tDep->SetLineColor(4);
    gr_DDVCSPlusBH_tDep->SetLineStyle(1);
    gr_DDVCSPlusBH_tDep->SetLineWidth(4);
    gr_DDVCSPlusBH_tDep->Draw("AL");
    gr_DDVCSPlusBH_tDep->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dtdx_{B}d#Phi dQ^{'2}[nb/GeV^{6}]");

    c1->Clear();
    TLegend *leg_DDVCS = new TLegend(0.6, 0.75, 0.92, 0.89);
    leg_DDVCS->SetBorderSize(0);
    leg_DDVCS->AddEntry(gr_DDVCSPlusBH_tDep, "DDVCS + BH");
    leg_DDVCS->AddEntry(gr_DDVCS_tDep, "DDVCS");
    TMultiGraph *mtgr_DDVCS = new TMultiGraph();
    mtgr_DDVCS->Add(gr_DDVCS_tDep);
    mtgr_DDVCS->Add(gr_DDVCSPlusBH_tDep);
    mtgr_DDVCS->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dtdx_{B}d#Phi dQ^{'2}[nb/GeV^{6}]");
    mtgr_DDVCS->Draw("AL");
    mtgr_DDVCS->GetXaxis()->SetLimits(0., 1.2);
    mtgr_DDVCS->SetMinimum(2.e-10);
    leg_DDVCS->Draw();
    c1->Print("Figs/DDVCS_10p6GeV_t_dep.pdf");
    c1->Print("Figs/DDVCS_10p6GeV_t_dep.png");
    c1->Print("Figs/DDVCS_10p6GeV_t_dep.root");
    
    TGraph *gr_DVCS_10p6 = new TGraph("DVCS_10p6GeV_xsec.dat", "%lg %*lg %lg");
    gr_DVCS_10p6->SetLineColor(2);
    gr_DVCS_10p6->SetLineStyle(8);
    gr_DVCS_10p6->SetLineWidth(4);
    
    TGraph *gr_DVCSPlusBH_10p6 = new TGraph("DVCSPlusBH_10p6GeV_xsec.dat", "%lg %*lg %lg");
    gr_DVCSPlusBH_10p6->SetLineColor(2);
    gr_DVCSPlusBH_10p6->SetLineStyle(1);
    gr_DVCSPlusBH_10p6->SetLineWidth(4);
    
    c1->Clear();
    TLegend *leg_DVCS = new TLegend(0.6, 0.75, 0.92, 0.89);
    leg_DVCS->SetBorderSize(0);
    leg_DVCS->AddEntry(gr_DVCSPlusBH_10p6, "DVCS + BH");
    leg_DVCS->AddEntry(gr_DVCS_10p6, "DVCS");
    TMultiGraph *mtgr_DVCS1 = new TMultiGraph();
    mtgr_DVCS1->Add(gr_DVCS_10p6);
    mtgr_DVCS1->Add(gr_DVCSPlusBH_10p6);
    mtgr_DVCS1->SetTitle("; -t [GeV^{2}]; d#sigma/dQ^{2}dx_{B}dtd#Phi [nb/GeV^{4}]");
    mtgr_DVCS1->Draw(":AL");
    mtgr_DVCS1->GetXaxis()->SetLimits(0., 1.2);
    mtgr_DVCS1->SetMinimum(4.e-4);
    leg_DVCS->Draw();
    c1->Print(Form("Figs/DVCS_10pxGeV_t_dep.pdf"));
    c1->Print(Form("Figs/DVCS_10pxGeV_t_dep.png"));
    c1->Print(Form("Figs/DVCS_10pxGeV_t_dep.root"));
}