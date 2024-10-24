/* 
 * File:   UpdateConstants.cc
 * Author: rafopar
 *
 * Created on September 11, 2024, 6:51 PM
 */

#include <cstdlib>

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

void RemoveOutlier(TGraph *);

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
    int run = atoi(argv[1]);

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "", 950, 950);
    c1->SetTopMargin(0.04);
    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.13);

    std::map<int, double> m_Eb;

    m_Eb[1] = 10.6;
    m_Eb[2] = 22;

    double Eb = m_Eb[run];


    //    c1->SetTopMargin(0.02);
    //    c1->SetRightMargin(0.03);
    //    c1->SetLeftMargin(0.14);


    TFile file_in(Form("AnaJPsi_Run_%d.root", run), "Read");

    TF1 *f_Eloss = new TF1("f_Eloss", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 1.05 * Eb);
    f_Eloss->SetNpx(4500);

    TF1 *f_ThCorr = new TF1("f_ThCorr", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0])) + [5]/((x-[0])*(x-[0])*(x-[0]))", 0., 80.);
    f_ThCorr->SetNpx(4500);

    TF1 *f_PhiCorr = new TF1("f_PhiCorr", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 0.3 * Eb);
    f_PhiCorr->SetNpx(4500);

    TH2D *h_dPP_P_mup1 = (TH2D*) file_in.Get("h_dPP_P_mup1");
    h_dPP_P_mup1->SetTitle("; P_{Rec}(#mu^{+}) [GeV] ;  (P_{Rec} - P_{MC})/P_{MC}");
    SliceFit(h_dPP_P_mup1, 0.01, 0.7 * Eb, 35, f_Eloss, Form("dPP_P_mup1_Run_%d", run));

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
    SliceFit(h_dPP_P_mum1, 0.01, 0.7 * Eb, 32, f_Eloss, Form("dPP_P_mum1_Run_%d", run));

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
    SliceFit(h_DeltaTheta_Theta_mup1, 3, 39, 32, f_ThCorr, Form("DeltaTheta_Theta_mup_Run_%d", run));
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
    SliceFit(h_DeltaTheta_Theta_mum1, 4, 39, 32, f_ThCorr, Form("DeltaTheta_Theta_mum_Run_%d", run));
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
    SliceFit(h_DeltaPhi_Pt_mup1, 0.2, 0.18*Eb, 18, f_PhiCorr, Form("DeltaPhi_Pt_mup_Run_%d", run));
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
    SliceFit(h_DeltaPhi_Pt_mum1, 0.2, 0.18 * Eb, 18, f_PhiCorr, Form("DeltaPhi_Pt_mum_Run_%d", run));
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

    int gr_ind = 0;
    for (int i = 0; i < nBins; i++) {

        double x1 = min_X + i*delta_X;
        double x2 = min_X + (i + 1) * delta_X;
        int bin1 = inp_hist->GetXaxis()->FindBin(x1);
        int bin2 = inp_hist->GetXaxis()->FindBin(x2);

        TH1D *h_tmp = (TH1D*) inp_hist->ProjectionY(Form("h_tmp_%d", i), bin1, bin2);
        h_tmp->Draw();
        if( h_tmp->GetEntries() < 25 ){
            continue;
        }
        double rms = h_tmp->GetRMS();
        sp1->Search(h_tmp, 2 * rms, "nobackground", 0.2);

        inp_hist->SetAxisRange(x1, x2);
        double avg_x = inp_hist->GetMean(1);
        inp_hist->GetXaxis()->UnZoom();

        c_tmp->Print(Form("Figs/debugPlots_%s.pdf", keyWord.c_str()));
        double *xx = sp1->GetPositionX();
        double *yy = sp1->GetPositionY();
        int n_peaks = sp1->GetNPeaks();

        if (n_peaks < 1) {
            continue;
        }

        /*
         *  It is possible sometime to get more than one peak.
         *  In these cases we will chose the one with Maximum Y value
         */
        double max = yy[0];
        int ind_max = 0;

        for (int j = 1; j < n_peaks; j++) {
            if (yy[j] > max) {
                max = yy[j];
                ind_max = j;
            }
        }

        double peak_pos = xx[ind_max];

        c_tmp->Clear();
        gr.SetPoint(gr_ind, avg_x, peak_pos);
        gr_ind = gr_ind + 1;
    }

    RemoveOutlier(&gr);

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
void FormatHist(T h) {

    //TH2D *h = new TH2D("dsd", "", 100, 0, 100, 200, 0, 100);

    h->SetTitleSize(0.05, "X");
    h->SetTitleSize(0.05, "Y");
    h->SetLabelSize(0.05, "X");
    h->SetLabelSize(0.05, "Y");
    h->GetYaxis()->SetTitleOffset(1.35);
}

/*
 * This method requires points of the graph to be sorted in X axis, i.e. x[i+1] > x[i] for any i
 */
void RemoveOutlier(TGraph *gr) {

    int N = gr->GetN();

    double x_Prev, x_Next, x_Cur;
    double y_Prev, y_Next, y_Cur;

    double rms = 0;

    /*
     * Loop over the grap elements (excluding the 1st and the last element), and calculate
     * how far is the current point wrt interpolated values obtained from previous and next points
     * The calculated the RMS of those values
     */
    for (int i = 2; i < N; i++) {
        gr->GetPoint(i - 2, x_Prev, y_Prev);
        gr->GetPoint(i, x_Cur, y_Cur);
        gr->GetPoint(i - 1, x_Next, y_Next);

        double y_interpol = y_Prev + (x_Cur - x_Prev)*(y_Next - y_Prev) / (x_Next - x_Prev);

        // The distance between the line passing through previous and next points and the given point
        double d = TMath::Abs((y_Next - y_Prev) * x_Cur - (x_Next - x_Prev) * y_Cur + x_Next * y_Prev - y_Next * x_Prev) / sqrt((y_Next - y_Prev)*(y_Next - y_Prev) + (x_Next - x_Prev)*(x_Next - x_Prev));

        //rms = rms + (y_Cur - y_interpol)*(y_Cur - y_interpol);
        rms = rms + d*d;
    }

    rms = sqrt(rms / double(N - 2));

    /*
     * Loop again over the elements, and throw points that are at least 5 rms away
     */

    bool pointsRemoved = true;

    while (pointsRemoved) {

        pointsRemoved = false;
        int N = gr-> GetN();

        for (int i = 2; i < N; i++) {
            gr->GetPoint(i - 2, x_Prev, y_Prev);
            gr->GetPoint(i, x_Cur, y_Cur);
            gr->GetPoint(i - 1, x_Next, y_Next);

            double y_interpol = y_Prev + (x_Cur - x_Prev)*(y_Next - y_Prev) / (x_Next - x_Prev);
            double d = TMath::Abs((y_Next - y_Prev) * x_Cur - (x_Next - x_Prev) * y_Cur + x_Next * y_Prev - y_Next * x_Prev) / sqrt((y_Next - y_Prev)*(y_Next - y_Prev) + (x_Next - x_Prev)*(x_Next - x_Prev));
            cout << "RMS = " << rms << "   x = " << x_Cur << "   d = " << d << "    diff = " << d / rms << " [rms]" << endl;

            if (d > 6 * rms ) {
                cout << "Point is being removed" << endl;
                gr->RemovePoint(i);
                pointsRemoved = true;
            }
        }
    }
}