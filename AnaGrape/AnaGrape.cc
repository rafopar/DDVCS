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

bool emAcc(TLorentzVector&);
bool mumAcc(TLorentzVector&);
bool mupAcc(TLorentzVector&);
bool protAcc(TLorentzVector&);

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

    const double InstLumin = 1.e37;
    const double secondPerDay = 3600 * 24;
    const double nDays = 100;
    const double totLumi = InstLumin * secondPerDay*nDays;
    const double pbn = 1.e-36;



    TRandom *rand = new TRandom();

    int run = atoi(argv[1]);

    // Beloe L_em is the one in the acceptance, and will be determined event by event
    TLorentzVector L_targ, L_beam, L_em, L_mum, L_mup, L_prot, L_mummup, L_q;

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
    tr11->Add(Form("../Data/Run_%d/grape_Run%d.root/h11", run, run));
    tr11->SetBranchAddress("xsec", &xsec);
    tr11->GetEntry(0);
    cout << "crs = " << xsec[0] << endl;



    TChain *tr1 = new TChain();
    tr1->Add(Form("../Data/Run_%d/grape_Run%d.root/h1", run, run));

    tr1->SetBranchAddress("px", &px);
    tr1->SetBranchAddress("py", &py);
    tr1->SetBranchAddress("pz", &pz);
    tr1->SetBranchAddress("pe", &pe);
    tr1->SetBranchAddress("pm", &pm);
    tr1->SetBranchAddress("kf", &kf);
    tr1->SetBranchAddress("sta", &sta);
    tr1->SetBranchAddress("npy", &npy);

    TFile file_out(Form("AnaGrape_%d.root", run), "Recreate");
    TH1D h_Minv1("h_Minv1", "", 200, 0., 4.);
    TH2D h_Qp2_vs_Q2_1("h_Qp2_vs_Q2_1", "", 200, 0., 1., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_2("h_Qp2_vs_Q2_2", "", 200, 0., 10., 200, 0., 10);
    TH2D h_thP_Qp2_1("h_thP_Qp2_1", "", 200, 0., 10., 200, 0., 90.);
    TH2D h_thP_Qp2_2("h_thP_Qp2_2", "", 200, 0., 10., 200, 0., 90.);
    TH2D h_thP_tM_1("h_thP_tM_1", "", 200, 0., 3.5, 200, 0., 90.);
    TH2D h_thP_tM_2("h_thP_tM_2", "", 200, 0., 3.5, 200, 0., 90.);
    TH2D h_thP_tM_3("h_thP_tM_3", "", 200, 0., 3.5, 200, 0., 90.);
    TH2D h_th_P_p1("h_th_P_p1", "", 200, 0., 4.5, 200, 0., 80);
    TH2D h_th_P_p2("h_th_P_p2", "", 200, 0., 4.5, 200, 0., 80);
    TH2D h_th_P_p3("h_th_P_p3", "", 200, 0., 4.5, 200, 0., 80);

    int iFile = 0;

    int nev = tr1->GetEntries();

    for (int i_ev = 0; i_ev < nev; i_ev++) {
        //for (int i_ev = 0; i_ev < 100; i_ev++) {
        if (i_ev % 50000 == 0) {
            cout.flush() << "Processed " << i_ev << "\r";
            //            if (i_ev % 5000000 == 0) {
            //                gDirectory->Write("hist_Dir", TObject::kOverwrite);
            //            }
        }

        tr1->GetEntry(i_ev);
        L_prot.SetPxPyPzE(px[10], py[10], pz[10], pe[10]);
        L_prot.RotateY(PI);
        L_em.SetPxPyPzE(px[11], py[11], pz[11], pe[11]);
        L_em.RotateY(PI);
        L_mup.SetPxPyPzE(px[12], py[12], pz[12], pe[12]);
        L_mup.RotateY(PI);
        L_mum.SetPxPyPzE(px[13], py[13], pz[13], pe[13]);
        L_mum.RotateY(PI);
        L_beam.SetPxPyPzE(px[1], py[1], pz[1], pe[1]);
        L_beam.RotateY(PI);
        L_targ.SetPxPyPzE(px[0], py[0], pz[0], pe[0]);
        L_targ.RotateY(PI);

        L_mummup = L_mum + L_mup;
        h_Minv1.Fill(L_mummup.M());

        L_q = L_em - L_beam;
        double tM = 2 * Mp * (L_prot.E() - Mp);
        double Qp2 = L_mummup.M2();
        double Q2 = -L_q.M2();

        double th_prot = L_prot.Theta() * radian;
        double p_prot = L_prot.P();

        h_Qp2_vs_Q2_1.Fill(Q2, Qp2);
        bool em_Acc = emAcc(L_em);
        bool mum_Acc = mumAcc(L_mum);
        bool mup_Acc = mupAcc(L_mup);
        bool prot_Acc = protAcc(L_mup);

        h_thP_Qp2_1.Fill(Qp2, th_prot);
        h_thP_tM_1.Fill(tM, th_prot);
        h_th_P_p1.Fill(p_prot, th_prot);

        if (em_Acc && mum_Acc && mup_Acc) {
            h_Qp2_vs_Q2_2.Fill(Q2, Qp2);
            h_thP_Qp2_2.Fill(Qp2, th_prot);
            h_thP_tM_2.Fill(tM, th_prot);
            h_th_P_p2.Fill(p_prot, th_prot);

            if (Qp2 > 2) {
                h_thP_tM_3.Fill(tM, th_prot);
                h_th_P_p3.Fill(p_prot, th_prot);
            }
        }
    }

    double weight = xsec[0] * pbn * totLumi / nev;

    TList *l_objs = gDirectory->GetList();

    for (TIter curObj = l_objs->begin(); curObj != l_objs->end(); curObj.Next()) {
        ((TH1*) * curObj)->Scale(weight);
    }

    gDirectory->Write();
    return 0;
}

bool emAcc(TLorentzVector& L) {

    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max;
}

bool mumAcc(TLorentzVector& L) {

    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.5; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max;
}

bool mupAcc(TLorentzVector& L) {

    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.5; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max;
}

bool protAcc(TLorentzVector& L) {

    const double th_max = 120; // deg
    const double th_min = 40.; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.5; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max;
}