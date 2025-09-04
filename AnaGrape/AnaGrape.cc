/* 
 * File:   ConvertLund.cc
 * Author: rafopar
 *
 * Created on March 29, 2021, 4:10 PM
 */

#include <fstream>
#include <iostream>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <cstdlib>

using namespace std;

bool emAcc(TLorentzVector&);
bool mumAcc(TLorentzVector&);
bool mupAcc(TLorentzVector&);
bool protAcc(TLorentzVector&);

TF1 *f_mumAccThP;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 2) {

        cout << "Please provide the 'Run' number " << endl;
        return 1;
    }
    f_mumAccThP = new TF1("f_mumAccThP", "[0] + [1]/(x-[2])", 0., 25);
    f_mumAccThP->SetParameters(5.07868, 18.1913, -0.120759);

    const double PI = TMath::Pi();
    const double r2d = TMath::RadToDeg();
    const double Mp = 0.9383;
    const double me = 0.00051;

    const double InstLumin = 1.e37;
    const double secondPerDay = 3600 * 24;
    const double nDays = 200;
    const double totLumi = InstLumin * secondPerDay*nDays;
    const double pbn = 1.e-36;
    const double Qp2Cut_EPIC = 1.2; // GeV

    const int nQp2bins = 4;
    const int nQ2bins = 4;
    const int nQp2bins_bin3 = 3; // The number of Qp2 bins in Kinematic bin3
    const double Qp2_edges[nQp2bins + 1] = {0.8, 1.6, 2.4, 3.2, 4.};
    const double Q2_edges[nQ2bins + 1] = {1., 1.6, 2.4, 3.2, 4.};
    const double Qp2_edges_bin2[nQp2bins + 1] = {1., 2., 3., 4., 6.};
    const double Qp2_edges_bin3[nQp2bins_bin3 + 1] = {1., 2., 3., 5.};

    // When the Qp2 is in between Qp2_Min and Qp2_Max, we will scan over Q2, i.e. will have nQ2bins, and for each fill phi_LH and "\xi vs x" histograms
    const double Qp2_Max = 3.;
    const double Qp2_Min = 2.5;
    const double Q2_Max = 2;
    const double Q2_Min = 1.5;
    const double Qp2_ResFreeCut = 4.;

    const double Q2_edges_new[nQ2bins + 1] = {1., 1.3, 1.7, 2.5, 4.};
    const double Qp2_edges_new[nQ2bins + 1] = {2., 2.3, 2.9, 3.7, 5.};

    TH1D h_Q2_Edges_new("h_Q2_Edges_new", "", nQ2bins, Q2_edges_new);
    TH1D h_Qp2_Edges_new("h_Qp2_Edges_new", "", nQp2bins, Qp2_edges_new);

    TH1D h_Qp2_Edges("h_Qp2_Edges", "", nQp2bins, Qp2_edges);
    TH1D h_Qp2_Edges_bin2("h_Qp2_Edges_bin2", "", nQp2bins, Qp2_edges_bin2);
    TH1D h_Qp2_Edges_bin3("h_Qp2_Edges_bin3", "", nQp2bins_bin3, Qp2_edges_bin3);
    TH1D h_Q2_Edges("h_Q2_Edges", "", nQ2bins, Q2_edges);



    TRandom *rand = new TRandom();

    int run = atoi(argv[1]);

    // Below L_em is the one in the acceptance, and will be determined event by event
    TLorentzVector L_targ, L_beam, L_em, L_mum, L_mup, L_prot, L_mummup, L_q, L_W;

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
    TH1D h_Minv2("h_Minv2", "", 200, 0., 4.);
    TH2D h_Qp2_vs_Q2_1("h_Qp2_vs_Q2_1", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_2("h_Qp2_vs_Q2_2", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_3("h_Qp2_vs_Q2_3", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_4("h_Qp2_vs_Q2_4", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_5("h_Qp2_vs_Q2_5", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_6("h_Qp2_vs_Q2_6", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_7("h_Qp2_vs_Q2_7", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_8("h_Qp2_vs_Q2_8", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_9("h_Qp2_vs_Q2_9", "", 200, 0., 10., 200, 0., 10);
    TH2D h_thP_Qp2_1("h_thP_Qp2_1", "", 200, 0., 10., 200, 0., 90.);
    TH2D h_thP_Qp2_2("h_thP_Qp2_2", "", 200, 0., 10., 200, 0., 90.);
    TH2D h_thP_tM_1("h_thP_tM_1", "", 200, 0., 3.5, 200, 0., 90.);
    TH2D h_thP_tM_2("h_thP_tM_2", "", 200, 0., 3.5, 200, 0., 90.);
    TH2D h_thP_tM_3("h_thP_tM_3", "", 200, 0., 3.5, 200, 0., 90.);
    TH2D h_th_P_p1("h_th_P_p1", "", 200, 0., 4.5, 200, 0., 80);
    TH2D h_th_P_p2("h_th_P_p2", "", 200, 0., 4.5, 200, 0., 80);
    TH2D h_th_P_p3("h_th_P_p3", "", 200, 0., 4.5, 200, 0., 80);

    TH2D h_th_P_em1("h_th_P_em1", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mum1("h_th_P_mum1", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mup1("h_th_P_mup1", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_em2("h_th_P_em2", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mum2("h_th_P_mum2", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mup2("h_th_P_mup2", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_em3("h_th_P_em3", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mum3("h_th_P_mum3", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mup3("h_th_P_mup3", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_Q2_xB1("h_Q2_xB1", "", 200, 0., 1., 200, 0., 10.);
    TH2D h_Q2_xB2("h_Q2_xB2", "", 200, 0., 1., 200, 0., 10.);
    TH2D h_Q2_xB3("h_Q2_xB3", "", 200, 0., 1., 200, 0., 10.);
    TH2D h_Q2_xB4("h_Q2_xB4", "", 200, 0., 1., 200, 0., 10.);

    TH2D h_xiprime_vs_xi1("h_xiprime_vs_xi1", "", 200, 0., 0.5, 200, -0.5, 0.5);
    TH2D h_Qp2_tM1("h_Qp2_tM1", "", 200, 0., 2., 200, 0., 8.);
    TH2D h_Qp2_tM2("h_Qp2_tM2", "", 200, 0., 2., 200, 0., 8.);
    TH2D h_Qp2_tM3("h_Qp2_tM3", "", 200, 0., 2., 200, 0., 8.);
    TH2D h_Qp2_tM4("h_Qp2_tM4", "", 200, 0., 2., 200, 0., 8.);

    TH2D h_xB_tM1("h_xB_tM1", "", 200, 0., 1.5, 200, 0., 1.);
    TH2D h_xB_tM2("h_xB_tM2", "", 200, 0., 1.5, 200, 0., 1.);
    TH2D h_xB_tM6("h_xB_tM6", "", 200, 0., 1.5, 200, 0., 1.);
    TH2D h_xB_tM7("h_xB_tM7", "", 200, 0., 1.5, 200, 0., 1.);
    TH2D h_xB_tM8("h_xB_tM8", "", 200, 0., 1.5, 200, 0., 1.);
    TH2D h_xB_tM9("h_xB_tM9", "", 200, 0., 1.5, 200, 0., 1.);

    TH2D h_xi_xxGPD1("h_xi_xxGPD1", "", 200, -1., 1., 200, 0., 1.);
    TH2D h_xi_xxGPD2("h_xi_xxGPD2", "", 200, -1., 1., 200, 0., 1.);
    TH2D h_xi_xxGPD7("h_xi_xxGPD7", "", 200, -1., 1., 200, 0., 1.);

    TH2D h_xi_xxGPD2_[nQp2bins];
    TH2D h_xi_xxGPD3_[nQp2bins];
    TH2D h_xi_xxGPD_bin2_[nQp2bins]; // xi_xxGPD in the kine bin 2
    TH2D h_xi_xxGPD_bin3_[nQp2bins_bin3]; // xi_xxGPD in the kine bin 3
    TH1D h_Phi_LH_bin1_Qp2_[nQp2bins];
    TH1D h_Phi_LH_bin1_Q2_[nQp2bins];
    TH1D h_Phi_LH_bin2_Qp2_[nQp2bins];
    TH1D h_Phi_LH_bin3_Qp2_[nQp2bins_bin3];

    for (int i = 0; i < nQp2bins; i++) {
        h_xi_xxGPD2_[i] = TH2D(Form("h_xi_xxGPD2_%d", i), "", 200, -0.5, 0.5, 200, 0., 0.5);
        h_xi_xxGPD3_[i] = TH2D(Form("h_xi_xxGPD3_%d", i), "", 200, -0.5, 0.5, 200, 0., 0.5);
        h_xi_xxGPD_bin2_[i] = TH2D(Form("h_xi_xxGPD_bin2_%d", i), "", 200, -0.5, 0.5, 200, 0., 0.5);

        h_Phi_LH_bin1_Qp2_[i] = TH1D(Form("h_Phi_LH_bin1_Qp2_%d", i), "", 12, 0, 360);
        h_Phi_LH_bin1_Q2_[i] = TH1D(Form("h_Phi_LH_bin1_Q2_%d", i), "", 12, 0, 360);
        h_Phi_LH_bin2_Qp2_[i] = TH1D(Form("h_Phi_LH_bin2_Qp2_%d", i), "", 12, 0, 360);
    }
    for (int i = 0; i < nQp2bins_bin3; i++) {
        h_xi_xxGPD_bin3_[i] = TH2D(Form("h_xi_xxGPD_bin3_%d", i), "", 200, -0.5, 0.5, 200, 0., 0.5);
        h_Phi_LH_bin3_Qp2_[i] = TH1D(Form("h_Phi_LH_bin3_Qp2_%d", i), "", 12, 0, 360);
    }


    TH2D h_xi_xxGPD_newQ2Scan_[nQ2bins]; // Scan over Q2 for a fixed Qp2
    TH1D h_Phi_LH_newQ2Scan_[nQ2bins]; // Scan over Q2 for a fixed Qp2
    TH2D h_Qp2_Q2_newQ2Scan_[nQ2bins]; // We need this to find the avg value for Michel's code
    TH2D h_tM_xB_newQ2Scan_[nQ2bins]; // We need this to find the avg value for Michel's code

    for (int i = 0; i < nQ2bins; i++) {
        h_xi_xxGPD_newQ2Scan_[i] = TH2D(Form("h_xi_xxGPD_newQ2Scan_%d", i), "", 200, -0.5, 0.5, 200, 0., 0.5);
        h_Phi_LH_newQ2Scan_[i] = TH1D(Form("h_Phi_LH_newQ2Scan_%d", i), "", 12, 0, 360);
        h_Qp2_Q2_newQ2Scan_[i] = TH2D(Form("h_Qp2_Q2_newQ2Scan_%d", i), "", 200, 0., 10., 200, 0., 10.);
        h_tM_xB_newQ2Scan_[i] = TH2D(Form("h_tM_xB_newQ2Scan_%d", i), "", 200, 0., 1.2, 200, 0., 1.);
    }

    TH2D h_xi_xxGPD_newQp2Scan_[nQp2bins]; // Scan over Qp2 for a fixed Q2
    TH1D h_Phi_LH_newQp2Scan_[nQp2bins]; // Scan over Qp2 for a fixed Q2
    TH2D h_Qp2_Q2_newQp2Scan_[nQ2bins]; // We need this to find the avg value for Michel's code
    TH2D h_tM_xB_newQp2Scan_[nQ2bins]; // We need this to find the avg value for Michel's code
    for (int i = 0; i < nQp2bins; i++) {
        h_xi_xxGPD_newQp2Scan_[i] = TH2D(Form("h_xi_xxGPD_newQp2Scan_%d", i), "", 200, -0.5, 0.5, 200, 0., 0.5);
        h_Phi_LH_newQp2Scan_[i] = TH1D(Form("h_Phi_LH_newQp2Scan_%d", i), "", 12, 0, 360);
        h_Qp2_Q2_newQp2Scan_[i] = TH2D(Form("h_Qp2_Q2_newQp2Scan_%d", i), "", 200, 0., 10., 200, 0., 10.);
        h_tM_xB_newQp2Scan_[i] = TH2D(Form("h_tM_xB_newQp2Scan_%d", i), "", 200, 0., 1.2, 200, 0., 1.);
    }

    // New Binning on \xi and x instead of Q2, Qp2, -t and xB

    const double bin0_x_max = -0.05;
    const double bin0_x_min = -0.09;
    const double bin0_xi_max = 0.245;
    const double bin0_xi_min = 0.225;

    const double bin1_x_max = 0.09;
    const double bin1_x_min = 0.05;
    const double bin1_xi_max = 0.245;
    const double bin1_xi_min = 0.225;

    const double t_Max_xi_xiStudy = 0.5;
    const int n_xi_x_bins = 2;
    const int n_xi_x_Q2bins = 2;


    TH2D h_xi_xxGPD_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH1D h_Phi_LH_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D h_Qp2_Q2_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D h_tM_xB_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 

    for (int i = 0; i < n_xi_x_bins; i++) {
        for (int j = 0; j < n_xi_x_bins; j++) {
            h_xi_xxGPD_xi_x_[i][j] = TH2D(Form("h_xi_xxGPD_xi_x_%d_%d", i, j), "", 200, -0.5, 0.5, 200, 0., 0.5);
            h_Phi_LH_xi_x_[i][j] = TH1D(Form("h_Phi_LH_xi_x_%d_%d", i, j), "", 12, 0, 360);
            h_Qp2_Q2_xi_x_[i][j] = TH2D(Form("h_Qp2_Q2_xi_x_%d_%d", i, j), "", 200, 0., 10., 200, 0., 10.);
            h_tM_xB_xi_x_[i][j] = TH2D(Form("h_tM_xB_xi_x_%d_%d", i, j), "", 200, 0., 1.2, 200, 0., 1.);
        }
    }



    int iFile = 0;

    int nev = tr1->GetEntries();
    //nev = 1500000;

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

        L_mummup = L_mum + L_mup; // The LorentzVector for the muon pair
        h_Minv1.Fill(L_mummup.M());

        L_q = L_em - L_beam; // The LorentzVector of Spacelaike photon
        TVector3 v3_beam_Eprime = L_beam.Vect().Cross(L_em.Vect());
        TVector3 v3_q_qprime = L_q.Vect().Cross(L_mummup.Vect());
        double Phi_LH = L_em.Vect().Dot(v3_q_qprime) > 0 ? v3_beam_Eprime.Angle(v3_q_qprime) * r2d : 360. - v3_beam_Eprime.Angle(v3_q_qprime) * r2d;
        double tM = 2 * Mp * (L_prot.E() - Mp); // Calculating the Mandelstam t
        double Qp2 = L_mummup.M2(); // Calculating the Q2 Prime as the invariant mass square of muon pair
        double Q2 = -L_q.M2(); // Calculating the spacelike Q2
        double nue = -L_q.E(); // New, defined as the energy of the spacelike photon
        double xB = Q2 / (2 * Mp * nue); // Calculating the xB
        double xiPrime = xB / (2 - xB); // Calculating the xiPrime as xB/(2-xB)
        double xi = xiPrime * (Q2 + Qp2) / (Q2); // Calculating xi as xiPrime*(Q2 + Qp2)/Q2


        //cout<<"xx_GPD = "<<xx_GPD<<endl;

        double xi_1Prime = (Q2 - Qp2 - tM / 2.) / (2 * Q2 / xB - Q2 - Qp2 - tM);
        double xi_1 = (Q2 + Qp2) / (2 * Q2 / xB - Q2 - Qp2 - tM);

        double xx_GPD = 2 * xiPrime - xi; // This is x that GPDs depend on defined as 2*xiPrime  - xi

        //cout<<xi<<"   "<<xi_1<<endl;
        //cout<<xiPrime<<"   "<<xi_1Prime<<endl;

        double th_em = L_em.Theta() * r2d;
        double p_em = L_em.P();

        double th_mup = L_mup.Theta() * r2d;
        double p_mup = L_mup.P();
        double th_mum = L_mum.Theta() * r2d;
        double p_mum = L_mum.P();

        double th_prot = L_prot.Theta() * r2d;
        double p_prot = L_prot.P();

        h_Qp2_vs_Q2_1.Fill(Q2, Qp2);
        bool em_Acc = emAcc(L_em);
        bool mum_Acc = mumAcc(L_mum);
        bool mup_Acc = mupAcc(L_mup);
        bool prot_Acc = protAcc(L_prot);

        h_thP_Qp2_1.Fill(Qp2, th_prot);
        h_thP_tM_1.Fill(tM, th_prot);
        h_th_P_p1.Fill(p_prot, th_prot);
        h_th_P_em1.Fill(p_em, th_em);
        h_th_P_mup1.Fill(p_mup, th_mup);
        h_th_P_mum1.Fill(p_mum, th_mum);
        h_xB_tM1.Fill(tM, xB);

        if (Q2 > Qp2) {
            h_xi_xxGPD1.Fill(xx_GPD, xi);
        }

        h_xiprime_vs_xi1.Fill(xi, xiPrime);
        h_Q2_xB1.Fill(xB, Q2);

        int Qp2Bin = h_Qp2_Edges.FindBin(Qp2) - 1;
        int Q2Bin = h_Q2_Edges.FindBin(Q2) - 1;

        h_Qp2_tM1.Fill(tM, Qp2);

        if (em_Acc && mum_Acc && mup_Acc /* && prot_Acc */) {

            h_Minv2.Fill(L_mummup.M());

            h_Qp2_vs_Q2_2.Fill(Q2, Qp2);
            h_thP_Qp2_2.Fill(Qp2, th_prot);
            h_thP_tM_2.Fill(tM, th_prot);
            h_xi_xxGPD2.Fill(xx_GPD, xi);
            h_th_P_p2.Fill(p_prot, th_prot);
            h_th_P_em2.Fill(p_em, th_em);
            h_th_P_mup2.Fill(p_mup, th_mup);
            h_th_P_mum2.Fill(p_mum, th_mum);
            h_Q2_xB2.Fill(xB, Q2);
            h_xB_tM2.Fill(tM, xB);

            if (Qp2 > Qp2_ResFreeCut) {
                h_Q2_xB4.Fill(xB, Q2);
            }

            if (Qp2 > Qp2Cut_EPIC) {
                h_th_P_p3.Fill(p_prot, th_prot);
                h_th_P_em3.Fill(p_em, th_em);
                h_th_P_mup3.Fill(p_mup, th_mup);
                h_th_P_mum3.Fill(p_mum, th_mum);
            }


            if (Qp2 > 0.3 && Qp2 < 0.6) {
                h_Q2_xB3.Fill(xB, Q2);
            }

            if (tM < 0.4) {
                h_Qp2_vs_Q2_6.Fill(Q2, Qp2);

                if (Qp2 > Qp2_Min && Qp2 < Qp2_Max) {
                    //cout<<"Qp2_Min = "<<Qp2_Min<<"    Qp2_Max = "<<Qp2_Max << "   Qp2 = " << Qp2 << "   Q2 = " << Q2 << endl;
                    int new_Q2bin = h_Q2_Edges_new.FindBin(Q2) - 1;

                    if (new_Q2bin >= 0 && new_Q2bin < nQ2bins) {

                        h_Phi_LH_newQ2Scan_[new_Q2bin].Fill(Phi_LH);
                        h_xi_xxGPD_newQ2Scan_[new_Q2bin].Fill(xx_GPD, xi);

                        h_Qp2_Q2_newQ2Scan_[new_Q2bin].Fill(Q2, Qp2);
                        h_tM_xB_newQ2Scan_[new_Q2bin].Fill(tM, xB);

                    }
                }

                if (Q2 > Q2_Min && Q2 < Q2_Max) {
                    int new_Qp2bin = h_Qp2_Edges_new.FindBin(Qp2) - 1;

                    if (new_Qp2bin >= 0 && new_Qp2bin < nQp2bins) {
                        h_Phi_LH_newQp2Scan_[new_Qp2bin].Fill(Phi_LH);
                        h_xi_xxGPD_newQp2Scan_[new_Qp2bin].Fill(xx_GPD, xi);
                        h_Qp2_Q2_newQp2Scan_[new_Qp2bin].Fill(Q2, Qp2);
                        h_tM_xB_newQp2Scan_[new_Qp2bin].Fill(tM, xB);
                    }
                }

            }


            if (tM < t_Max_xi_xiStudy) {

                int bin_xi_x = -1;
                int bin_Q2 = -1;

                if (xx_GPD > bin0_x_min && xx_GPD < bin0_x_max && xi > bin0_xi_min && xi < bin0_xi_max) {
                    bin_xi_x = 0;

                    bin_Q2 = Qp2 > 4.5 - Q2 / 1.65 ? 1 : 0;

                } else if (xx_GPD > bin1_x_min && xx_GPD < bin1_x_max && xi > bin1_xi_min && xi < bin1_xi_max) {
                    bin_xi_x = 1;

                    bin_Q2 = Qp2 > 8.2 - Q2 / 0.55 ? 1 : 0;
                }



                if (bin_xi_x >= 0 && bin_Q2 >= 0) {
                    h_Phi_LH_xi_x_[bin_xi_x][bin_Q2].Fill(Phi_LH);
                    h_xi_xxGPD_xi_x_[bin_xi_x][bin_Q2].Fill(xx_GPD, xi);
                    h_Qp2_Q2_xi_x_[bin_xi_x][bin_Q2].Fill(Q2, Qp2);
                    h_tM_xB_xi_x_[bin_xi_x][bin_Q2].Fill(tM, xB);
                }

            }


            //h_xi_xxGPD8
            if (tM < 0.4) {
                h_Qp2_vs_Q2_7.Fill(Q2, Qp2);

                if (Q2 > 1 && Q2 < 2.05 && Qp2 > 3.5 && Qp2 < 4. && xB < 0.07) {
                    h_xi_xxGPD7.Fill(xx_GPD, xi);
                    h_xB_tM7.Fill(tM, xB);
                }
            }

            if (xx_GPD > -0.08 && xx_GPD < -0.06 && xi > 0.23 && xi < 0.24) {
                h_Qp2_vs_Q2_8.Fill(Q2, Qp2);
                h_xB_tM8.Fill(tM, xB);
            } else if (xx_GPD > 0.06 && xx_GPD < 0.08 && xi > 0.23 && xi < 0.24) {
                h_Qp2_vs_Q2_9.Fill(Q2, Qp2);
                h_xB_tM9.Fill(tM, xB);
            }


            if (xB > 0.12 && xB < 0.22 && tM > 0.1 && tM < 0.4) {
                h_Qp2_tM2.Fill(tM, Qp2);
                h_Qp2_vs_Q2_3.Fill(Q2, Qp2);

                if (Q2 > 2. && Q2 < 3.) {
                    if (Qp2Bin >= 0 && Qp2Bin < nQp2bins) {
                        h_xi_xxGPD2_[Qp2Bin].Fill(xx_GPD, xi);
                        h_Phi_LH_bin1_Qp2_[Qp2Bin].Fill(Phi_LH);
                    }
                }

                if (Qp2 > 2 && Qp2 < 3.) {
                    if (Q2Bin >= 0 && Q2Bin < nQ2bins) {
                        h_xi_xxGPD3_[Q2Bin].Fill(xx_GPD, xi);
                        h_Phi_LH_bin1_Q2_[Q2Bin].Fill(Phi_LH);
                    }
                }
            }


            // ================= Bin 2 ======================
            if (xB > 0.01 && xB < 0.05) {
                h_Qp2_tM3.Fill(tM, Qp2);

                if (tM > 0.05 && tM < 0.4) {
                    h_Qp2_vs_Q2_4.Fill(Q2, Qp2);

                    int Qp2bn_bin2 = h_Qp2_Edges_bin2.FindBin(Qp2) - 1;
                    if (Qp2bn_bin2 >= 0 && Qp2bn_bin2 < nQp2bins) {
                        h_xi_xxGPD_bin2_[Qp2bn_bin2].Fill(xx_GPD, xi);
                        h_Phi_LH_bin2_Qp2_[Qp2bn_bin2].Fill(Phi_LH);
                    }
                }

            }

            if (xB > 0.05 && xB < 0.12) {
                h_Qp2_tM4.Fill(tM, Qp2);
                if (tM > 0.05 && tM < 0.4) {
                    h_Qp2_vs_Q2_5.Fill(Q2, Qp2);

                    int Qp2bn_bin3 = h_Qp2_Edges_bin3.FindBin(Qp2) - 1;
                    if (Qp2bn_bin3 >= 0 && Qp2bn_bin3 < nQp2bins_bin3) {
                        h_xi_xxGPD_bin3_[Qp2bn_bin3].Fill(xx_GPD, xi);
                        h_Phi_LH_bin3_Qp2_[Qp2bn_bin3].Fill(Phi_LH);
                    }
                }
            }



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

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max &&
            L.Theta() * TMath::RadToDeg() > f_mumAccThP->Eval(L.P());
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
    const double P_min = 0.3; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max;
}