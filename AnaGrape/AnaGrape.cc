/* 
 * File:   ConvertLund.cc
 * Author: rafopar
 *
 * Created on March 29, 2021, 4:10 PM
 */

#include <fstream>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TRandom.h>
#include <TLorentzVector.h>

#include <cstdlib>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 2) {

        cout << "Please provide the 'Run' number " << endl;
        return 1;
    }

    const double PI = TMath::Pi();
    const double radian = TMath::RadToDeg();
    const double Mp = 0.9383;
    const double me = 0.00051;
    const double vz_max = 2.5; // cm
    const double vz_min = -2.5; // cm

    const int n_perfile = 10000; // # of events per Lund file

    TRandom *rand = new TRandom();

    int run = atoi(argv[1]);

    // Beloe L_em is the one in the acceptance, and will be determined event by event
    TLorentzVector L_targ, L_beam, L_em, L_em1, L_ep, L_em2, L_prot, L_emep;

    Double_t xsec[2];

    Float_t px[20];
    Float_t py[20];
    Float_t pz[20];
    Float_t pe[20];
    Float_t pm[20];
    Float_t kf[20];
    Float_t sta[20];
    Float_t npy[20];

    TChain *tr11 = new TChain();
    tr11->Add(Form("Grape_Run%d/grape_*.root/h11", run));
    tr11->SetBranchAddress("xsec", &xsec);
    tr11->GetEntry(0);
    cout << "crs = " << xsec[0] << endl;

    TChain *tr1 = new TChain();
    tr1->Add(Form("Grape_Run%d/grape_*.root/h1", run));

    tr1->SetBranchAddress("px", &px);
    tr1->SetBranchAddress("py", &py);
    tr1->SetBranchAddress("pz", &pz);
    tr1->SetBranchAddress("pe", &pe);
    tr1->SetBranchAddress("pm", &pm);
    tr1->SetBranchAddress("kf", &kf);
    tr1->SetBranchAddress("sta", &sta);
    tr1->SetBranchAddress("npy", &npy);

    int iFile = 0;
    ofstream lund_Out(Form("../Grape_Lund/RunCard_%d/grape_Run_%d_file_%d", run, run, iFile));

    int nev = tr1->GetEntries();

    for (int i_ev = 0; i_ev < nev; i_ev++) {
        //for (int i_ev = 0; i_ev < 100; i_ev++) {
        if (i_ev % 50000 == 0) {
            cout.flush() << "Processed " << i_ev << "\r";
            //            if (i_ev % 5000000 == 0) {
            //                gDirectory->Write("hist_Dir", TObject::kOverwrite);
            //            }
        }
            if ((i_ev + 1) % n_perfile == 0) {

                lund_Out.close();
                iFile++;
                //Lund_out.open(Form("JPsi_gen_%d.txt", file_number), ofstream::out);
                lund_Out.open(Form("../Grape_Lund/RunCard_%d/grape_Run_%d_file_%d", run, run, iFile));

            }

        tr1->GetEntry(i_ev);
        L_prot.SetPxPyPzE(px[10], py[10], pz[10], pe[10]);
        L_prot.RotateY(PI);
        L_em1.SetPxPyPzE(px[11], py[11], pz[11], pe[11]);
        L_em1.RotateY(PI);
        L_ep.SetPxPyPzE(px[12], py[12], pz[12], pe[12]);
        L_ep.RotateY(PI);
        L_em2.SetPxPyPzE(px[13], py[13], pz[13], pe[13]);
        L_em2.RotateY(PI);
        L_beam.SetPxPyPzE(px[1], py[1], pz[1], pe[1]);
        L_beam.RotateY(PI);
        L_targ.SetPxPyPzE(px[0], py[0], pz[0], pe[0]);
        L_targ.RotateY(PI);

        // we will generate a random vz along the target.
        double vz = rand->Uniform(vz_min, vz_max);


        lund_Out << 4 << setw(5) << 1 << setw(5) << 1 << setw(5) << 0 << setw(9) << 0.85 << setw(5) << 11 << setw(10) << L_beam.E() << setw(9) << 2212 << setw(5) << 0 << setw(5) << 0 << endl;
        lund_Out << 1 << setw(5) << 0 << setw(5) << 1 << setw(9) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << L_prot.Px() << setw(15) << L_prot.Py() << setw(15) << L_prot.Pz() << setw(15) << L_prot.E() << setw(15) << L_prot.M() << setw(5) << 0 << setw(5) << setw(5) << 0. << setw(15) << vz << endl;
        lund_Out << 2 << setw(5) << 0 << setw(5) << 1 << setw(9) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << L_em1.Px() << setw(15) << L_em1.Py() << setw(15) << L_em1.Pz() << setw(15) << L_em1.E() << setw(15) << L_em1.M() << setw(5) << 0 << setw(5) << 0. << setw(15) << vz << endl;
        lund_Out << 3 << setw(5) << 0 << setw(5) << 1 << setw(9) << -11 << setw(5) << 0 << setw(5) << 0 << setw(15) << L_ep.Px() << setw(15) << L_ep.Py() << setw(15) << L_ep.Pz() << setw(15) << L_ep.E() << setw(15) << L_ep.M() << setw(5) << 0 << setw(5) << 0. << setw(15) << vz << endl;
        lund_Out << 4 << setw(5) << 0 << setw(5) << 1 << setw(9) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << L_em2.Px() << setw(15) << L_em2.Py() << setw(15) << L_em2.Pz() << setw(15) << L_em2.E() << setw(15) << L_em2.M() << setw(5) << 0 << setw(5) << 0. << setw(15) << vz << endl;
    }

    return 0;
}

