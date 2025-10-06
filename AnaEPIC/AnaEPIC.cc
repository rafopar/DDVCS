/* 
 * File:   AnaEPIC.cc
 * Author: rafopar
 *
 * Created on September 18, 2023, 11:26 AM
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

    int run = -1;
    if (argc != 3) {

        cout << "Please provide the Generator type (grape or epic) and 'Run' number " << endl;
        cout << "As an example ./AnaEPIC grape 9" << endl;
        cout << "exiting" << endl;
        return 1;
    }
    f_mumAccThP = new TF1("f_mumAccThP", "[0] + [1]/(x-[2])", 0., 25);
    f_mumAccThP->SetParameters(5.07868, 18.1913, -0.120759);

    std::string genType = argv[1];
    run = atoi(argv[2]);

    cout << "Generator is " << genType << endl;

    std::map<int, double> m_epicXsec;
    std::map<int, double> m_grapeXsec;
    const double InstLumin = 1.e37;
    const double secondPerDay = 3600 * 24;
    const double nDays = 100;
    const double totLumi = InstLumin * secondPerDay*nDays;
    const double nb = 1.e-33;
    const double pb = 1.e-36;

    const double PI = TMath::Pi();
    const double r2d = TMath::RadToDeg();
    const double Eb = 10.2;
    const double Mp = 0.9383;

    const double MinvCut = 1.21;
    
    int nPart;
    const int nMaxPart = 250;
    int index[nMaxPart];
    int pid[nMaxPart];
    int type[nMaxPart];
    int parentInd[nMaxPart];
    int daughtInd[nMaxPart];
    double t_live[nMaxPart];
    double px[nMaxPart];
    double py[nMaxPart];
    double pz[nMaxPart];
    double E[nMaxPart];
    double m[nMaxPart];
    double vx[nMaxPart];
    double vy[nMaxPart];
    double vz[nMaxPart];

    m_epicXsec[10] = 0.0077862966745603 * nb; // nb
    m_epicXsec[2] = 0.00074946419241075 * nb; // nb
    m_epicXsec[3] = 0.00193917558968409 * nb; // nb
    m_epicXsec[4] = 0.00193238282227725 * nb; // nb
    m_grapeXsec[9] = 0.263586 * pb; // pb

    // Below L_em is the one in the acceptance, and will be determined event by event
    TLorentzVector L_targ, L_beam, L_em, L_mum, L_mup, L_prot, L_mummup, L_q, L_W;

    L_targ.SetPxPyPzE(0., 0., 0., Mp);
    L_beam.SetPxPyPzE(0., 0., Eb, Eb);


    TTree *tr1 = new TTree();
    TFile *file_In = new TFile(Form("../Data/%s/%s_Run_%d.root", genType.c_str(), genType.c_str(), run));
    tr1 = (TTree*) file_In->Get("tr1");

    tr1->SetBranchAddress("nPart", &nPart);
    tr1->SetBranchAddress("index", index);
    tr1->SetBranchAddress("t_live", t_live);
    tr1->SetBranchAddress("type", type);
    tr1->SetBranchAddress("pid", pid);
    tr1->SetBranchAddress("parentInd", parentInd);
    tr1->SetBranchAddress("daughtInd", daughtInd);
    tr1->SetBranchAddress("px", &px[0]);
    tr1->SetBranchAddress("py", py);
    tr1->SetBranchAddress("pz", pz);
    tr1->SetBranchAddress("E", E);
    tr1->SetBranchAddress("m", m);
    tr1->SetBranchAddress("vx", vx);
    tr1->SetBranchAddress("vy", vy);
    tr1->SetBranchAddress("vz", vz);

    TFile *file_out = new TFile(Form("Ana_%s_%d.root", genType.c_str(), run), "Recreate");
    TH1D h_Minv1("h_Minv1", "", 200, 0., 4.);
    TH1D h_Minv2("h_Minv2", "", 200, 0., 4.);
    TH2D h_Qp2_vs_Q2_1("h_Qp2_vs_Q2_1", "", 200, 0., 10., 200, 0., 10);
    TH2D h_Qp2_vs_Q2_2("h_Qp2_vs_Q2_2", "", 200, 0., 10., 200, 0., 10);
    TH2D h_th_P_em2("h_th_P_em2", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mum2("h_th_P_mum2", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_mup2("h_th_P_mup2", "", 200, 0., 11., 200, 0., 40.);
    TH2D h_th_P_p2("h_th_P_p2", "", 200, 0., 4.5, 200, 0., 80);
    TH2D h_Q2_xB1("h_Q2_xB1", "", 200, 0., 1.05, 200, 0., 10.);
    TH2D h_Q2_xB2("h_Q2_xB2", "", 200, 0., 1.05, 200, 0., 10.);
    TH1D h_y1("h_y1", "", 200, 0., 1.);
    TH1D h_y2("h_y2", "", 200, 0., 1.);
    TH1D h_tM1("h_tM1", "", 200, 0., 9);
    TH1D h_tM2("h_tM2", "", 200, 0., 9);


    int nEv = tr1->GetEntries();
    //  nEv = 10;

    for (int iev = 0; iev < nEv; iev++) {

        if (iev % 50000 == 0) {
            cout.flush() << "Processed " << iev << "\r";
        }

        tr1->GetEntry(iev);

        int ind_em = -1;
        int ind_prot = -1;
        int ind_mum = -1;
        int ind_mup = -1;


        if (nPart < 4) {
            cout << "npart =" << nPart << " this should not happen!!!" << endl;
        }

        for (int ipart = 0; ipart < nPart; ipart++) {
            switch (pid[ipart]) {
                case 2212:
                    ind_prot = ipart;
                    break;
                case 11:
                    ind_em = ipart;
                    break;
                case 13:
                    ind_mum = ipart;
                    break;
                case -13:
                    ind_mup = ipart;
                    break;
                default:
                    break;
            }

        }



        L_em.SetPxPyPzE(px[ind_em], py[ind_em], pz[ind_em], E[ind_em]);
        L_prot.SetPxPyPzE(px[ind_prot], py[ind_prot], pz[ind_prot], E[ind_prot]);
        L_mum.SetPxPyPzE(px[ind_mum], py[ind_mum], pz[ind_mum], E[ind_mum]);
        L_mup.SetPxPyPzE(px[ind_mup], py[ind_mup], pz[ind_mup], E[ind_mup]);

        //cout<<ind_em<<"   "<<ind_prot<<"   "<<ind_mum<<"   "<<ind_mup<<endl;

        L_mummup = L_mum + L_mup;
        double minv = L_mummup.M();
        h_Minv1.Fill(minv);

        L_q = L_em - L_beam; // The LorentzVector of Spacelaike photon
        TVector3 v3_beam_Eprime = L_beam.Vect().Cross(L_em.Vect());
        TVector3 v3_q_qprime = L_q.Vect().Cross(L_mummup.Vect());
        double Phi_LH = L_em.Vect().Dot(v3_q_qprime) > 0 ? v3_beam_Eprime.Angle(v3_q_qprime) * r2d : 360. - v3_beam_Eprime.Angle(v3_q_qprime) * r2d;
        double tM = 2 * Mp * (L_prot.E() - Mp); // Calculating the Mandelstam t
        double Qp2 = L_mummup.M2(); // Calculating the Q2 Prime as the invariant mass square of muon pair
        double Q2 = -L_q.M2(); // Calculating the spacelike Q2
        double nue = -L_q.E(); // New, defined as the energy of the spacelike photon
        double y = nue / Eb;
        double xB = Q2 / (2 * Mp * nue); // Calculating the xB
        double xiPrime = xB / (2 - xB); // Calculating the xiPrime as xB/(2-xB)
        double xi = xiPrime * (Q2 + Qp2) / (Q2); // Calculating xi as xiPrime*(Q2 + Qp2)/Q2

        // ====== Pawel wants me to compare th_CM for grape and EPIC
        
        
        h_Qp2_vs_Q2_1.Fill(Q2, Qp2);
        h_Q2_xB1.Fill(xB, Q2);
        h_y1.Fill(y);

        h_tM1.Fill(tM);

        bool em_Acc = emAcc(L_em);
        bool mum_Acc = mumAcc(L_mum);
        bool mup_Acc = mupAcc(L_mup);
        bool prot_Acc = protAcc(L_prot);

        if (em_Acc && mum_Acc && mup_Acc && minv > MinvCut ) {
            h_Minv2.Fill(minv);
            h_Qp2_vs_Q2_2.Fill(Q2, Qp2);
            h_Q2_xB2.Fill(xB, Q2);
            h_y2.Fill(y);
            h_tM2.Fill(tM);

            h_th_P_em2.Fill(L_em.P(), L_em.Theta() * r2d);
            h_th_P_mum2.Fill(L_mum.P(), L_mum.Theta() * r2d);
            h_th_P_mup2.Fill(L_mup.P(), L_mup.Theta() * r2d);
            h_th_P_p2.Fill(L_prot.P(), L_prot.Theta() * r2d);

        }

    }

    TList *l_objs = gDirectory->GetList();

    double weight = 0;

    cout << "Kuku 0 Weight = " << weight << endl;
    if (genType.compare("grape") == 0) {
        weight = m_grapeXsec[run] * totLumi / double(nEv);
        cout << "Kuku 1 Weight = " << weight << endl;
    } else if (genType.compare("epic") == 0) {
        weight = m_epicXsec[run] * totLumi / double(nEv);
        cout << "Kuku 2 Weight = " << weight << endl;
    }

    cout << "Kuku 3 Weight = " << weight << endl;
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
