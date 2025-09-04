/* 
 * File:   UpdateConstants.cc
 * Author: rafopar
 *
 * Created on September 11, 2024, 6:51 PM
 */

#include <cstdlib>

#include <map>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <fstream>

using namespace std;

void SliceFit(TH2D* inp_hist, double min_X, double max_X, int nBins, TF1 *f_Fit, std::string keyWord);

void DumpFitParameters(TF1*, std::string keyWord);
template <typename T> void FormatHist(T h);

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Please provde the run number. Exiting..." << endl;
        exit(1);
    }

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.04);
    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.13);

    //    c1->SetTopMargin(0.02);
    //    c1->SetRightMargin(0.03);
    //    c1->SetLeftMargin(0.14);
    
    

    int run = atoi(argv[1]);

    std::map<int, double> m_Eb;
    
    m_Eb[17] = 10.6;
    m_Eb[32] = 11.;
    m_Eb[32001] = 11.;
    m_Eb[32002] = 11.;
    m_Eb[32003] = 11.;
    m_Eb[32004] = 11.;
    m_Eb[32005] = 11.;
    m_Eb[18] = 22;
    
    double Eb = m_Eb[run];
    
    TFile file_in(Form("AnaDDVCS_Run_%d.root", run), "Read");

    TF1 *f_Eloss = new TF1("f_Eloss", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 80.);
    f_Eloss->SetNpx(4500);

    TF1 *f_ThCorr = new TF1("f_ThCorr", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0])) + [5]/((x-[0])*(x-[0])*(x-[0]))", 0., 80.);
    f_ThCorr->SetNpx(4500);

    TF1 *f_PhiCorr = new TF1("f_PhiCorr", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 80.);
    f_PhiCorr->SetNpx(4500);
    
    TH2D *h_dPP_P_mup1 = (TH2D*) file_in.Get("h_dPP_P_mup1");
    h_dPP_P_mup1->SetTitle("; P_{Rec}(#mu^{+}) [GeV] ;  (P_{Rec} - P_{MC})/P_{MC}");
    SliceFit(h_dPP_P_mup1, 0.2, 0.6*Eb, 30, f_Eloss, Form("dPP_P_mup1_Run_%d", run));

    c1->cd();
    FormatHist(h_dPP_P_mup1);

    h_dPP_P_mup1->Draw();
    f_Eloss->Draw("Same");
    c1->Print(Form("Figs/mup_ElossFit_Run_%d.pdf", run));
    c1->Print(Form("Figs/mup_ElossFit_Run_%d.png", run));
    c1->Print(Form("Figs/mup_ElossFit_Run_%d.root", run));
    //DumpFitParameters(f_Eloss, Form("mup_ElossFunc_Pars_Run_%d", run));

    TH2D *h_dPP_P_mum1 = (TH2D*) file_in.Get("h_dPP_P_mum1");
    h_dPP_P_mum1->SetTitle("; P_{Rec}(#mu^{-}) [GeV] ;  (P_{Rec} - P_{MC})/P_{MC}");
    SliceFit(h_dPP_P_mum1, 0.06, 0.55*Eb, 30, f_Eloss, Form("dPP_P_mum1_Run_%d", run));

    c1->cd();
    FormatHist(h_dPP_P_mum1);
    h_dPP_P_mum1->Draw();
    f_Eloss->Draw("Same");
    c1->Print(Form("Figs/mum_ElossFit_Run_%d.pdf", run));
    c1->Print(Form("Figs/mum_ElossFit_Run_%d.png", run));
    c1->Print(Form("Figs/mum_ElossFit_Run_%d.root", run));
    //DumpFitParameters(f_Eloss, Form("mum_ElossFunc_Pars_Run_%d", run));

    TH2D *h_DeltaTheta_Theta_mup1 = (TH2D*) file_in.Get("h_DeltaTheta_Theta_mup1");
    h_DeltaTheta_Theta_mup1->SetTitle("; #theta_{Rec} [deg]; #theta_{Rec} - #theta_{MC} [deg]");
    SliceFit(h_DeltaTheta_Theta_mup1, 3, 40, 18, f_ThCorr, Form("DeltaTheta_Theta_mup_Run_%d", run));
    c1->cd();
    h_DeltaTheta_Theta_mup1->Draw();
    FormatHist(h_DeltaTheta_Theta_mup1);
    f_ThCorr->Draw("Same");
    c1->Print(Form("Figs/DeltaTheta_Theta_mup1_Run_%d.pdf", run));
    c1->Print(Form("Figs/DeltaTheta_Theta_mup1_Run_%d.png", run));
    c1->Print(Form("Figs/DeltaTheta_Theta_mup1_Run_%d.root", run));
    //DumpFitParameters(f_ThCorr, Form("mup_DeltaTheta_Corr_Run_%d", run));

    TH2D *h_DeltaTheta_Theta_mum1 = (TH2D*) file_in.Get("h_DeltaTheta_Theta_mum1");
    h_DeltaTheta_Theta_mum1->SetTitle("; #theta_{Rec} [deg]; #theta_{Rec} - #theta_{MC} [deg]");
    SliceFit(h_DeltaTheta_Theta_mum1, 6, 40, 18, f_ThCorr, Form("DeltaTheta_Theta_mum_Run_%d", run));
    c1->cd();
    h_DeltaTheta_Theta_mum1->Draw();
    FormatHist(h_DeltaTheta_Theta_mum1);
    f_ThCorr->Draw("Same");
    c1->Print(Form("Figs/DeltaTheta_Theta_mum1_Run_%d.pdf", run));
    c1->Print(Form("Figs/DeltaTheta_Theta_mum1_Run_%d.png", run));
    c1->Print(Form("Figs/DeltaTheta_Theta_mum1_Run_%d.root", run));
    //DumpFitParameters(f_ThCorr, Form("mum_DeltaTheta_Corr_Run_%d", run));

    TH2D *h_DeltaPhi_Pt_mup1 = (TH2D*) file_in.Get("h_DeltaPhi_Pt_mup1");
    h_DeltaPhi_Pt_mup1->SetTitle("; P_{t} [GeV]; #phi_{Rec} - #phi_{MC} [deg]");
    SliceFit(h_DeltaPhi_Pt_mup1, 0, 0.15*Eb, 14, f_PhiCorr, Form("DeltaPhi_Pt_mup_Run_%d", run));
    c1->cd();
    h_DeltaPhi_Pt_mup1->Draw();
    FormatHist(h_DeltaPhi_Pt_mup1);
    f_PhiCorr->Draw("Same");
    c1->Print(Form("Figs/DeltaPhi_Pt_mup1_Run_%d.pdf", run));
    c1->Print(Form("Figs/DeltaPhi_Pt_mup1_Run_%d.png", run));
    c1->Print(Form("Figs/DeltaPhi_Pt_mup1_Run_%d.root", run));
    //DumpFitParameters(f_PhiCorr, Form("mup_DeltaPhi_Corr_Run_%d", run));

    TH2D *h_DeltaPhi_Pt_mum1 = (TH2D*) file_in.Get("h_DeltaPhi_Pt_mum1");
    h_DeltaPhi_Pt_mum1->SetTitle("; P_{t} [GeV]; #phi_{Rec} - #phi_{MC} [deg]");
    SliceFit(h_DeltaPhi_Pt_mum1, 0, 0.15*Eb, 14, f_PhiCorr, Form("DeltaPhi_Pt_mum_Run_%d", run));
    c1->cd();
    h_DeltaPhi_Pt_mum1->Draw();
    FormatHist(h_DeltaPhi_Pt_mum1);
    f_PhiCorr->Draw("Same");
    c1->Print(Form("Figs/DeltaPhi_Pt_mum1_Run_%d.pdf", run));
    c1->Print(Form("Figs/DeltaPhi_Pt_mum1_Run_%d.png", run));
    c1->Print(Form("Figs/DeltaPhi_Pt_mum1_Run_%d.root", run));
    //DumpFitParameters(f_PhiCorr, Form("mum_DeltaPhi_Corr_Run_%d", run));

    return 0;
}

void SliceFit(TH2D* inp_hist, double min_X, double max_X, int nBins, TF1 *f_Fit, std::string keyWord) {

    TSpectrum *sp1 = new TSpectrum();

    TGraph gr;
    gr.SetMarkerColor(2);
    gr.SetMarkerStyle(20);

    double delta_X = (max_X - min_X) / double(nBins);

    TCanvas *c_tmp = new TCanvas("c_tmp", "", 950, 950);

    c_tmp->Print(Form("Figs/debugPlots_%s.pdf[", keyWord.c_str()));

    for (int i = 0; i < nBins; i++) {

        double x1 = min_X + i*delta_X;
        double x2 = min_X + (i + 1) * delta_X;
        int bin1 = inp_hist->GetXaxis()->FindBin(x1);
        int bin2 = inp_hist->GetXaxis()->FindBin(x2);

        TH1D *h_tmp = (TH1D*) inp_hist->ProjectionY(Form("h_tmp_%d", i), bin1, bin2);
        h_tmp->Draw();
        double rms = h_tmp->GetRMS();
        sp1->Search(h_tmp, 1 * rms, "nobackground", 0.2);

        inp_hist->SetAxisRange(x1, x2);
        double avg_x = inp_hist->GetMean(1);
        inp_hist->GetXaxis()->UnZoom();

        c_tmp->Print(Form("Figs/debugPlots_%s.pdf", keyWord.c_str()));
        double *xx = sp1->GetPositionX();
        double *yy = sp1->GetPositionY();
        int n_peaks = sp1->GetNPeaks();

        /*
         *  It is possible sometime to get more than one peak.
         *  In these cases we will chose the one with Maximum Y value
         */
        double max = yy[0];
        int ind_max = 0;

        for (int i = 1; i < n_peaks; i++) {
            if (yy[i] > max) {
                max = yy[i];
                ind_max = i;
            }
        }

        double peak_pos = xx[ind_max];

        c_tmp->Clear();
        gr.SetPoint(i, avg_x, peak_pos);
    }

    gr.Fit(f_Fit, "MeV", "", min_X, max_X);
    inp_hist->Draw();
    gr.Draw("P Same");
    c_tmp->Print(Form("Figs/debugPlots_%s.pdf", keyWord.c_str()));

    c_tmp->Print(Form("Figs/debugPlots_%s.pdf]", keyWord.c_str()));
}

void DumpFitParameters(TF1* f, std::string keyWord) {
    ofstream out_dat(Form("Pars/%s.dat", keyWord.c_str()), std::ios::out | std::ios::trunc);

    int nPar = f->GetNpar();

    for (int i = 0; i < nPar; i++) {
        out_dat << i << "  " << f->GetParameter(i) << endl;
    }

    out_dat.close();
}

template <typename T>
void FormatHist(T h){

    //TH2D *h = new TH2D("dsd", "", 100, 0, 100, 200, 0, 100);
    
    h->SetTitleSize(0.05, "X");
    h->SetTitleSize(0.05, "Y");
    h->SetLabelSize(0.05, "X");
    h->SetLabelSize(0.05, "Y");
    h->GetYaxis()->SetTitleOffset(1.35);
}
