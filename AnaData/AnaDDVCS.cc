/* 
 * File:   AnaDDVCS.cc
 * Author: rafopar
 *
 * Created on August 17, 2023, 8:40 PM
 */

#include <cstdlib>
#include <iostream>

#include <TF1.h>
#include <TKey.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TChain.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <ElectronDDVCSKine.h>
#include <RecParticle.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
#include <vector>

using namespace std;

bool emAcc(TLorentzVector &);

bool mumAcc(TLorentzVector &);

bool mupAcc(TLorentzVector &);

bool protAcc(TLorentzVector &);

TF1 *f_mumAccThP;

/*
 *
 */
int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Please provide the run number." << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    int run = atoi(argv[1]);
    char inputFile[256];

    std::map<int, double> m_TotXsec;
    std::map<int, double> m_Eb;
    std::map<int, double> m_totCharge;
    std::map<int, double> m_BeamPol;
    std::map<int, bool> m_isMC;

    sprintf(inputFile, "inpHipoFiles/ememep_Run_%d.hipo", run);

    /*
    * Let's read the RunSettings.dat file and initialize variables
    * that are specified in the RunSettings.dat
    */

    std::ifstream inp_config("RunSettings.dat");
    if (!inp_config.is_open()) {
        std::cerr << "Error opening the config file! \n Exiting..." << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(inp_config, line)) {
        size_t commentPos = line.find("#"); // Or line.find("#"); for hash comments

        // If a comment is found, truncate the string at that position
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string entry;
        std::vector<std::string> v_entries;

        while (iss >> entry) {
            v_entries.push_back(entry);
        }

        int arun = std::stoi(v_entries[0]);
        double axSec = std::stof(v_entries[1]);
        double aEb = std::stof(v_entries[2]);
        double aCarge = std::stof(v_entries[3]);
        double aBeamPol = std::stof(v_entries[4]);
        int aisMC = std::stoi(v_entries[5]);

        m_TotXsec[arun] = axSec;
        m_Eb[arun] = aEb;
        m_totCharge[arun] = aCarge;
        m_BeamPol[arun] = aBeamPol;
        m_isMC[arun] = aisMC;
    }



    /*
     * Check if the run exist in the RunSerrings.dat
     */

    if ( m_isMC.find(run) == m_isMC.end() ) {
        cerr<<"The run "<<run<<" is not in the RunSettings.dat. Exiting..."<<endl;
        exit(1);
    }


    bool isMC = m_isMC[run];
    double x_sec = m_TotXsec[run];
    double Eb = m_Eb[run];
    double totCarge = m_totCharge[run];
    double beamPol = m_BeamPol[run];

    f_mumAccThP = new TF1("f_mumAccThP", "[0] + [1]/(x-[2])", 0., 25);
    f_mumAccThP->SetParameters(5.07868, 18.1913, -0.120759);

    const double Hist_Pmax = 11.;
    const int nMax_samePID = 20;
    const double Mp = 0.9383;
    const double chi2PID_Prot_Max = 10.;
    const double chi2PID_Prot_Min = -10.;
    const int HTCC_TYPE = 15;
    const int DET_FTOF = 12;
    const int layer1b = 2; // panel 1a = 1, panely2 = 3
    const double HTCC_MomThr = 4.8;
    const double PCalEmin = 0.07;
    const double PCalVmin = 10.;
    const double PCalWmin = 10.;
    const double vzMax = 1.;
    const double vzMin = -8.;
    const double v_dtCut = 1.5;
    const double Mmis_Max = 1.2;
    const double Mmis_Min = 0.8;

    const double P_prot_Max = 1.5; // MC shows that the majority of protons have momentum less than 0.8 GeV
    const double dP_P_cut = 0.3;

    const double Mx2ReactionCut = 0.1; // GeV^2
    const double MinvMinCut = 0.9; // To consider for DDVCS analysis, the Minv should be above this cut.

    const double PmisCut = 1.4; // In GeV. Momenta of about 90% of protons is below this value
    const double thMinCut = 15.; // Majority of protons fly above this value

    const double light_Speed = 29.9792458;

    /*
     * Defining constants for normalizing against S19 Lumi
     */
    const double rho = 0.07085; // gm*cm^{-3}
    const double l = 5.; // cm
    const double N_A = 6.022e23;
    const double em_charge = 1.6e-19; //
    const double N_targ = rho * l * N_A;
    constexpr double mili = 0.001;
    constexpr double inv_pb = 1e-36;
    double N_e_onTarg = (totCarge*mili)/em_charge;
    double Lumi = N_targ*N_e_onTarg*inv_pb;

    const double MxRecoil_Max = 1.1; // GeV
    const double MxRecoil_Min = 0.88; // GeV


    TF1 *f_Separation = new TF1("f_Separation", "pol6", 0., 9);
    f_Separation->SetParameters(24.7129, -19.3813, 11.4942, -3.88646, 0.738393, -0.0723433, 0.00279077);

    ElectronDDVCSKine ddvcs_kine;

    ddvcs_kine.SetEb(Eb);
    ddvcs_kine.SetMtarg(Mp);

    ElectronDDVCSKine ddvcs_kine_Noprot;
    ddvcs_kine_Noprot.SetEb(Eb);
    ddvcs_kine_Noprot.SetMtarg(Mp);

    TLorentzVector L_beam;
    L_beam.SetPxPyPzE(0., 0., Eb, Eb);
    TLorentzVector L_targ;
    L_targ.SetPxPyPzE(0., 0., 0., Mp);

    TObjArray Hlist(0);

    TFile *file_out = new TFile(Form("HistData/Hists_DDVCS_Run_%d.root", run), "Recreate");
    TH1D h_nphe_em1("h_nphe_em1", "", 200, 0., 35.);
    TH1D h_nphe_ep1("h_nphe_ep1", "", 200, 0., 35.);

    TH1D h_vz_em1("h_vz_em1", "", 200, -25., 35.);
    TH1D h_vz_ep1("h_vz_ep1", "", 200, -25., 35.);
    TH1D h_vz_prot1("h_vz_prot1", "", 200, -25., 35.);
    TH1D h_chi2PID_em1("h_chi2PID_em1", "", 400, -200, 200.);
    TH1D h_chi2PID_ep1("h_chi2PID_ep1", "", 400, -200, 200.);
    TH1D h_chi2PID_prot1("h_chi2PID_prot1", "", 400, -200, 200.);

    TH2D h_th_P_em1("h_th_P_em1", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_ep1("h_th_P_ep1", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_prot1("h_th_P_prot1", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);

    TH2D h_th_P_em1_2("h_th_P_em1_2", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em2_2("h_th_P_em2_2", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em2("h_th_P_em2", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_ep2("h_th_P_ep2", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_prot2("h_th_P_prot2", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);

    TH2D h_th_P_em1_3("h_th_P_em1_3", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em2_3("h_th_P_em2_3", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em3("h_th_P_em3", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_ep3("h_th_P_ep3", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_prot3("h_th_P_prot3", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);

    TH2D h_th_P_em1_4("h_th_P_em1_4", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em2_4("h_th_P_em2_4", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em4("h_th_P_em4", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_ep4("h_th_P_ep4", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);

    TH2D h_th_P_em1_5("h_th_P_em1_5", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em2_5("h_th_P_em2_5", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_em5("h_th_P_em5", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);
    TH2D h_th_P_ep5("h_th_P_ep5", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 55.);

    TH2D h_beta_P_prot1("h_beta_P_prot1", "", 200, 0., 0.5 * Hist_Pmax, 200, 0., 1.5);
    TH2D h_beta_P_prot4("h_beta_P_prot4", "", 200, 0., 0.5 * Hist_Pmax, 200, 0., 1.5);
    TH1D h_E_PCal_em1("h_E_PCal_em1", "", 200, 0., 1.);
    TH1D h_E_PCal_ep1("h_E_PCal_ep1", "", 200, 0., 1.);

    TH1D h_lUPCal_em1("h_lUPCal_em1", "", 200, 0., 450.);
    TH1D h_lVPCal_em1("h_lVPCal_em1", "", 200, 0., 450.);
    TH1D h_lWPCal_em1("h_lWPCal_em1", "", 200, 0., 450.);
    TH1D h_lUPCal_ep1("h_lUPCal_ep1", "", 200, 0., 450.);
    TH1D h_lVPCal_ep1("h_lVPCal_ep1", "", 200, 0., 450.);
    TH1D h_lWPCal_ep1("h_lWPCal_ep1", "", 200, 0., 450.);
    TH2D h_SF_PCal_ECin_ep1("h_SF_PCal_ECin_ep1", "", 200, 0., 0.3, 200, 0., 0.3);
    TH2D h_SF_PCal_ECin_ep2("h_SF_PCal_ECin_ep2", "", 200, 0., 0.3, 200, 0., 0.3);
    TH2D h_n_emep1("h_n_emep1", "", 11, -0.5, 10.5, 11, -0.5, 10.5);
    TH1D h_n_prot1("h_n_prot1", "", 11, -0.5, 10.5);
    TH1D h_n_charged1("h_n_charged1", "", 11, -0.5, 10.5);

    TH1D h_Mmis1("h_Mmis1", "", 200, 0., 4.);
    TH2D h_th_VS_Mmis1("h_th_VS_Mmis1", "", 200, 0., 65., 200, 0., 4.);
    TH2D h_th_VS_Mmis2("h_th_VS_Mmis2", "", 200, 0., 65., 200, 0., 4.);
    TH2D h_th_VS_Mmis3("h_th_VS_Mmis3", "", 200, 0., 65., 200, 0., 4.);
    TH2D h_th_VS_Mmis4("h_th_VS_Mmis4", "", 200, 0., 65., 200, 0., 4.);

    TH2D h_Mmis_PMis1("h_Mmis_PMis1", "", 200, 0., 8., 200, 0., 4.);
    TH2D h_Mmis_PMis2("h_Mmis_PMis2", "", 200, 0., 8., 200, 0., 4.);
    TH2D h_Mmis_PMis3("h_Mmis_PMis3", "", 200, 0., 8., 200, 0., 4.);
    TH2D h_th_P_Mis1("h_th_P_Mis1", "", 200, 0., 8., 200, 0., 65);
    TH2D h_th_P_Mis2("h_th_P_Mis2", "", 200, 0., 8., 200, 0., 65);
    TH2D h_th_P_Mis3("h_th_P_Mis3", "", 200, 0., 8., 200, 0., 65);
    TH1D h_Mmis2("h_Mmis2", "", 200, 0., 4.);
    TH1D h_Mmis3("h_Mmis3", "", 200, 0., 4.);
    TH1D h_Mmis4("h_Mmis4", "", 200, 0., 4.);
    TH1D h_Mmis5("h_Mmis5", "", 200, 0., 4.);
    TH1D h_Mmis6("h_Mmis6", "", 200, 0., 4.);
    TH2D h_Minv12_1("h_Minv12_1", "", 200, 0., 2.5, 200, 0., 2.5);
    TH2D h_Minv12_2("h_Minv12_2", "", 200, 0., 2.5, 200, 0., 2.5);
    TH2D h_Minv12_3("h_Minv12_3", "", 200, 0., 2.5, 200, 0., 2.5);
    TH2D h_Minv12_4("h_Minv12_4", "", 200, 0., 2.5, 200, 0., 2.5);
    TH2D h_Minv12_5("h_Minv12_5", "", 200, 0., 2.5, 200, 0., 2.5);
    TH2D h_Minv12_NoProt3("h_Minv12_NoProt3", "", 200, 0., 2.5, 200, 0., 2.5);
    TH1D h_vt_Diff_em1("h_vt_Diff_em1", "", 200, -0.5, 0.5);
    TH2D h_vt_Diff_emep1("h_vt_Diff_emep1", "", 200, -2., 2., 200, -2., 2.);
    TH2D h_vz_Diff_emep1("h_vz_Diff_emep1", "", 200, -10., 10., 200, -10., 10.);
    TH1D h_vt_Diff_em2("h_vt_Diff_em2", "", 200, -18.5, 18.5);
    TH1D h_vt_Diff_em1_ep("h_vt_Diff_em1_ep", "", 200, -18.5, 18.5);
    TH1D h_vt_Diff_em2_ep("h_vt_Diff_em2_ep", "", 200, -18.5, 18.5);
    TH2D h_MC_Memep_12_1("h_MC_Memep_12_1", "", 200, 0., 3., 200, 0., 3.);

    TH2D h_Pmiss_Pprot4("h_Pmiss_Pprot4", "", 200, 0., 0.5*Eb, 200, 0., 0.5*Eb);
    TH2D h_dP_P_prot4("h_dP_P_prot4", "", 200, 0., 0.5*Eb, 200, -1, 1);

    TH1D h_Mmis4Part4("h_Mmis4Part4", "", 200, -1.5, 1.5);

    TH2D h_xi_xxGPD_1_NoProt3("h_xi_xxGPD_1_NoProt3", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_xi_xxGPD_2_NoProt3("h_xi_xxGPD_2_NoProt3", "", 200, -0.5, 0.5, 200, 0., 0.5);

    TH2D h_xi_xxGPD_1_4("h_xi_xxGPD_1_4", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_xi_xxGPD_2_4("h_xi_xxGPD_2_4", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_xi_xxGPD_1_5("h_xi_xxGPD_1_5", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_xi_xxGPD_2_5("h_xi_xxGPD_2_5", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_xi_xxGPD_1_6("h_xi_xxGPD_1_6", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_xi_xxGPD_2_6("h_xi_xxGPD_2_6", "", 200, -0.5, 0.5, 200, 0., 0.5);



    TH1D h_Phi_LH1_NoProt3("h_Phi_LH1_NoProt3", "", 12, 0., 360);
    TH1D h_Phi_LH2_NoProt3("h_Phi_LH2_NoProt3", "", 12, 0., 360);
    TH1D h_Phi_LH1_helWeighted_NoProt3("h_Phi_LH1_helWeighted_NoProt3", "", 12, 0., 360);
    TH1D h_Phi_LH2_helWeighted_NoProt3("h_Phi_LH2_helWeighted_NoProt3", "", 12, 0., 360);

    TH1D h_Phi_LH1_4("h_Phi_LH1_4", "", 12, 0., 360);
    TH1D h_Phi_LH2_4("h_Phi_LH2_4", "", 12, 0., 360);
    TH1D h_Phi_LH1_5("h_Phi_LH1_5", "", 12, 0., 360);
    TH1D h_Phi_LH2_5("h_Phi_LH2_5", "", 12, 0., 360);
    TH1D h_Phi_LH1_6("h_Phi_LH1_6", "", 12, 0., 360);
    TH1D h_Phi_LH2_6("h_Phi_LH2_6", "", 12, 0., 360);
    TH1D h_Phi_LH1_helWeighted_6("h_Phi_LH1_helWeighted_6", "", 12, 0., 360);
    TH1D h_Phi_LH2_helWeighted_6("h_Phi_LH2_helWeighted_6", "", 12, 0., 360);



    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;

    int evCounter = 0;

    hipo::bank bRunConf(factory.getSchema("RUN::config"));
    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));
    hipo::bank bRecSC(factory.getSchema("REC::Scintillator"));
    hipo::bank bRecEV(factory.getSchema("REC::Event"));
    hipo::bank bMCPart(factory.getSchema("MC::Particle"));

    int ind_em[nMax_samePID];
    int ind_ep[nMax_samePID];
    int ind_prot[nMax_samePID];
    int sec_em[nMax_samePID];
    int sec_ep[nMax_samePID];

    try {
        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            if (evCounter > 400000) {
                break;
            }

            if (evCounter % 10000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

            event.getStructure(bRecPart);
            event.getStructure(bRecCalo);
            event.getStructure(bRecCC);
            event.getStructure(bRecSC);
            event.getStructure(bRecEV);
            event.getStructure(bRunConf);

            map<int, int> ind_HTCC;
            map<int, int> ind_PCal;
            map<int, int> ind_ECin;
            map<int, int> ind_ECout;
            map<int, int> ind_FTOF;

            int helicity = bRecEV.getInt("helicity", 0);

            int nPart = bRecPart.getRows();
            for (int i_part = 0; i_part < nPart; i_part++) {
                // ==== Before assigning, index, all indexes are initialized to -1, this way we can check, whether
                // ==== that particular detector is related to the particle "i_part"
                ind_HTCC[i_part] = -1;
                ind_PCal[i_part] = -1;
                ind_ECin[i_part] = -1;
                ind_ECout[i_part] = -1;
                ind_FTOF[i_part] = -1;

                int nCherenkov = bRecCC.getRows();
                // =============================== HTCC ===============================
                for (int i_cc = 0; i_cc < nCherenkov; i_cc++) {
                    // Want only HTCC for now
                    if (bRecCC.getInt("detector", i_cc) == HTCC_TYPE) {
                        if (bRecCC.getInt("pindex", i_cc) == i_part) {
                            ind_HTCC[i_part] = i_cc;
                        }
                    }
                }

                int nTOF = bRecSC.getRows();

                // ===================== FTOF ========================
                for (int itof = 0; itof < nTOF; itof++) {
                    if (bRecSC.getInt("detector", itof) == DET_FTOF && bRecSC.getInt("layer", itof) == layer1b) {
                        if (bRecSC.getInt("pindex", itof) == i_part) {
                            ind_FTOF[i_part] = itof;
                        }
                    }
                }

                int nCal = bRecCalo.getRows();

                // =============================== PCal, ECin, ECout ===============================
                for (int i_cal = 0; i_cal < nCal; i_cal++) {
                    if (bRecCalo.getInt("pindex", i_cal) == i_part) {
                        int layer = bRecCalo.getInt("layer", i_cal);

                        if (layer == 1) {
                            ind_PCal[i_part] = i_cal;
                        } else if (layer == 4) {
                            ind_ECin[i_part] = i_cal;
                        } else if (layer == 7) {
                            ind_ECout[i_part] = i_cal;
                        }
                    }
                }
            }

            int n_em = 0;
            int n_ep = 0;
            int n_prot = 0;
            int n_charged = 0;


            for (int ipart = 0; ipart < nPart; ipart++) {
                RecParticle recp(bRecPart, bRecCalo, bRecCC, ipart, ind_PCal[ipart], ind_ECin[ipart], ind_ECout[ipart],
                                 ind_HTCC[ipart]);

                if (recp.pid() == 11 && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {
                    h_vz_em1.Fill(recp.vz());
                    h_chi2PID_em1.Fill(recp.chi2pid());
                    h_nphe_em1.Fill(recp.npheHTCC());

                    h_th_P_em1.Fill(recp.p(), recp.th());
                    h_E_PCal_em1.Fill(recp.energyPCal());

                    h_lUPCal_em1.Fill(recp.luPCal());
                    h_lVPCal_em1.Fill(recp.lvPCal());
                    h_lWPCal_em1.Fill(recp.lwPCal());

                    bool isPCalEmin = recp.energyPCal() > PCalEmin;
                    bool isPCalVmin = recp.lvPCal() > PCalVmin;
                    bool isPCalWmin = recp.lvPCal() > PCalWmin;
                    bool isVz = recp.vz() > vzMin && recp.vz() < vzMax;

                    if (isPCalEmin && isPCalVmin && isPCalWmin && isVz) {
                        ind_em[n_em] = ipart;
                        sec_em[n_em] = recp.phi() / 60;
                        n_em = n_em + 1;
                    }
                } else if (recp.pid() == -11 && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {
                    h_vz_ep1.Fill(recp.vz());
                    h_chi2PID_ep1.Fill(recp.chi2pid());
                    h_th_P_ep1.Fill(recp.p(), recp.th());
                    h_E_PCal_ep1.Fill(recp.energyPCal());

                    h_nphe_ep1.Fill(recp.npheHTCC());

                    h_lUPCal_ep1.Fill(recp.luPCal());
                    h_lVPCal_ep1.Fill(recp.lvPCal());
                    h_lWPCal_ep1.Fill(recp.lwPCal());

                    h_SF_PCal_ECin_ep1.Fill(recp.SFPCal(), recp.SFECin());

                    if (recp.p() > HTCC_MomThr) {
                        h_SF_PCal_ECin_ep2.Fill(recp.SFPCal(), recp.SFECin());
                    }

                    bool isPCalEmin = recp.energyPCal() > PCalEmin;
                    bool isPCalVmin = recp.lvPCal() > PCalVmin;
                    bool isPCalWmin = recp.lvPCal() > PCalWmin;

                    if (((recp.energyECin() < 0.001 && recp.SFPCal() > 0.11) || (
                             recp.energyECin() >= 0.001 && (recp.SFPCal() + recp.SFECin()) > 0.17)) && isPCalEmin &&
                        isPCalVmin && isPCalWmin) {
                        ind_ep[n_ep] = ipart;
                        sec_ep[n_ep] = recp.phi() / 60;
                        n_ep = n_ep + 1;
                    }
                } else if (recp.pid() == 2212) {
                    h_vz_prot1.Fill(recp.vz());
                    h_chi2PID_prot1.Fill(recp.chi2pid());
                    h_th_P_prot1.Fill(recp.p(), recp.th());
                    h_beta_P_prot1.Fill(recp.p(), recp.beta());

                    if (recp.chi2pid() > chi2PID_Prot_Min && recp.chi2pid() < chi2PID_Prot_Max) {
                        ind_prot[n_prot] = ipart;
                        n_prot = n_prot + 1;
                    }
                }
            }

            if (n_em == 2 && n_ep == 1) {
                TLorentzVector L_em1, L_em2, L_ep;

                RecParticle part_em1(bRecPart, bRecCalo, bRecCC, ind_em[0], ind_PCal[ind_em[0]], ind_ECin[ind_em[0]],
                                     ind_ECout[ind_em[0]], ind_HTCC[ind_em[0]]);

                RecParticle part_em2(bRecPart, bRecCalo, bRecCC, ind_em[1], ind_PCal[ind_em[1]], ind_ECin[ind_em[1]],
                                     ind_ECout[ind_em[1]], ind_HTCC[ind_em[1]]);
                RecParticle part_ep(bRecPart, bRecCalo, bRecCC, ind_ep[0], ind_PCal[ind_ep[0]], ind_ECin[ind_ep[0]],
                                    ind_ECout[ind_ep[0]], ind_HTCC[ind_ep[0]]);

                double v_t_em2 = bRecSC.getFloat("time", ind_FTOF[ind_em[1]]) - bRecSC.getFloat(
                                     "path", ind_FTOF[ind_em[1]]) / light_Speed;
                double v_t_ep = bRecSC.getFloat("time", ind_FTOF[ind_ep[0]]) - bRecSC.getFloat(
                                    "path", ind_FTOF[ind_ep[0]]) / light_Speed;


                if ( (part_em2.th() - f_Separation->Eval(part_em2.p()) ) > (part_em1.th() - f_Separation->Eval(part_em1.p()) ) ) {
                    L_em1.SetPxPyPzE(part_em1.px(), part_em1.py(), part_em1.pz(), part_em1.p());
                    L_em2.SetPxPyPzE(part_em2.px(), part_em2.py(), part_em2.pz(), part_em2.p());
                }else {
                    L_em1.SetPxPyPzE(part_em2.px(), part_em2.py(), part_em2.pz(), part_em2.p());
                    L_em2.SetPxPyPzE(part_em1.px(), part_em1.py(), part_em1.pz(), part_em1.p());
                }

                L_ep.SetPxPyPzE(part_ep.px(), part_ep.py(), part_ep.pz(), part_ep.p());

                ddvcs_kine_Noprot.SetKineem1em2ep(&L_em1, &L_em2, &L_ep);

                TLorentzVector L_mis = L_beam + L_targ - (L_em1 + L_em2 + L_ep);
                TLorentzVector L_emep1 = L_em2 + L_ep;
                TLorentzVector L_emep2 = L_em1 + L_ep;

                double m_emep1 = L_emep1.M();
                double m_emep2 = L_emep2.M();

                double th_mis = L_mis.Theta() * TMath::RadToDeg();
                double Mmis = L_mis.M();
                double Pmis = L_mis.P();

                double MxRecoil = ddvcs_kine_Noprot.GetMxRecoil();

                h_Mmis1.Fill(Mmis);
                h_Minv12_1.Fill(m_emep1, m_emep2);

                h_th_VS_Mmis1.Fill(th_mis, Mmis);
                h_Mmis_PMis1.Fill(Pmis, Mmis);
                h_th_P_Mis1.Fill(Pmis, th_mis);

                h_th_P_em1_2.Fill(part_em1.p(), part_em1.th());
                h_th_P_em2_2.Fill(part_em2.p(), part_em2.th());
                h_th_P_em2.Fill(part_em1.p(), part_em1.th());
                h_th_P_em2.Fill(part_em2.p(), part_em2.th());
                h_th_P_ep2.Fill(part_ep.p(), part_ep.th());


                if (Mmis > Mmis_Min && Mmis < Mmis_Max) {
                    h_Minv12_2.Fill(m_emep1, m_emep2);
                }

                h_vt_Diff_em1.Fill(part_em1.vt() - part_em2.vt());
                h_vt_Diff_em2.Fill(part_em1.vt() - v_t_em2);
                h_vt_Diff_emep1.Fill(part_em1.vt() - part_em2.vt(), part_em1.vt() - part_ep.vt());
                h_vz_Diff_emep1.Fill(part_em1.vz() - part_em2.vz(), part_em1.vz() - part_ep.vz());

                if (TMath::Abs(part_em1.vt() - v_t_em2) < v_dtCut) {
                    h_vt_Diff_em1_ep.Fill(part_em1.vt() - v_t_ep);
                    h_vt_Diff_em2_ep.Fill(v_t_em2 - v_t_ep);

                    if ( TMath::Abs(part_em1.vt() - v_t_ep) > v_dtCut || TMath::Abs(v_t_em2 - v_t_ep) > v_dtCut ) {
                        continue;
                    }

                    if (m_emep1 > 0.25 && m_emep2 > 0.25) {
                        // This is to skip low mass events which seems are not produced at the target
                        //if (m_emep1 > 1. && m_emep2 > 1.) { // This is to skip low mass events which seems are not produced at the target
                        h_Mmis3.Fill(Mmis);
                        h_th_VS_Mmis3.Fill(th_mis, Mmis);
                        h_Mmis_PMis3.Fill(Pmis, Mmis);
                        h_th_P_Mis3.Fill(Pmis, th_mis);

                        h_th_P_em1_3.Fill(part_em1.p(), part_em1.th());
                        h_th_P_em2_3.Fill(part_em2.p(), part_em2.th());
                        h_th_P_em3.Fill(part_em1.p(), part_em1.th());
                        h_th_P_em3.Fill(part_em2.p(), part_em2.th());
                        h_th_P_ep3.Fill(part_ep.p(), part_ep.th());

                        if ( MxRecoil > MxRecoil_Min && MxRecoil < MxRecoil_Max ) {

                            double xx_GPD1_NoProt = ddvcs_kine_Noprot.GetXX_GPD_1();
                            double xx_GPD2_NoProt = ddvcs_kine_Noprot.GetXX_GPD_2();
                            double xi1_NoProt = ddvcs_kine_Noprot.GetXi_1();
                            double xi2_NoProt = ddvcs_kine_Noprot.GetXi_2();
                            double phi_LH1_NoProt = ddvcs_kine_Noprot.GetPhi_LH_1();
                            double phi_LH2_NoProt = ddvcs_kine_Noprot.GetPhi_LH_2();

                            double Minv1_NoProt = ddvcs_kine_Noprot.GetMinv_1();
                            double Minv2_NoProt = ddvcs_kine_Noprot.GetMinv_2();

                            h_Minv12_NoProt3.Fill(Minv1_NoProt, Minv2_NoProt);

                            if ( Minv1_NoProt > MinvMinCut ) {
                                h_xi_xxGPD_1_NoProt3.Fill(xx_GPD1_NoProt, xi1_NoProt);
                                h_xi_xxGPD_2_NoProt3.Fill(xx_GPD2_NoProt, xi2_NoProt);
                                h_Phi_LH1_NoProt3.Fill(phi_LH1_NoProt);
                                h_Phi_LH1_helWeighted_NoProt3.Fill(phi_LH1_NoProt, helicity/beamPol);
                            }

                            if ( Minv2_NoProt > MinvMinCut ) {
                                h_Phi_LH2_NoProt3.Fill(phi_LH2_NoProt);
                                h_Phi_LH2_helWeighted_NoProt3.Fill(phi_LH2_NoProt, helicity/beamPol);
                            }
                        }

                        if (Pmis < PmisCut) {
                            h_th_VS_Mmis4.Fill(th_mis, Mmis);
                            h_th_P_em1_4.Fill(part_em1.p(), part_em1.th());
                            h_th_P_em2_4.Fill(part_em2.p(), part_em2.th());
                            h_th_P_em4.Fill(part_em1.p(), part_em1.th());
                            h_th_P_em4.Fill(part_em2.p(), part_em2.th());
                            h_th_P_ep4.Fill(part_ep.p(), part_ep.th());

                            if (th_mis > thMinCut && Mmis > Mmis_Min && Mmis < Mmis_Max) {
                                h_Minv12_3.Fill(m_emep1, m_emep2);
                                h_th_P_em1_5.Fill(part_em1.p(), part_em1.th());
                                h_th_P_em2_5.Fill(part_em2.p(), part_em2.th());
                                h_th_P_em5.Fill(part_em1.p(), part_em1.th());
                                h_th_P_em5.Fill(part_em2.p(), part_em2.th());
                                h_th_P_ep5.Fill(part_ep.p(), part_ep.th());
                            }
                        }

                        if (n_prot == 1) {
                            h_Mmis4.Fill(Mmis);
                            RecParticle part_p(bRecPart, bRecCalo, bRecCC, ind_prot[0], ind_PCal[ind_prot[0]], ind_ECin[ind_prot[0]], ind_ECout[ind_prot[0]], ind_HTCC[ind_prot[0]]);
                            h_beta_P_prot4.Fill(part_p.p(), part_p.beta());

                            TLorentzVector L_prot;
                            L_prot.SetPxPyPzE(part_p.px(), part_p.py(), part_p.pz(), sqrt( part_p.p()*part_p.p() + Mp*Mp) );
                            TLorentzVector L_mis4part = L_beam + L_targ - (L_em1 + L_em2 + L_ep + L_prot);

                            ddvcs_kine.SetKin4FSParticles(&L_em1, &L_em2, &L_ep, &L_prot);

                            double Mx2_4part = L_mis4part.M2();
                            h_Mmis4Part4.Fill(Mx2_4part);

                            double xi1 = ddvcs_kine.GetXi_1();
                            double xi2 = ddvcs_kine.GetXi_2();
                            double xx_GPD1 = ddvcs_kine.GetXX_GPD_1();
                            double xx_GPD2 = ddvcs_kine.GetXX_GPD_2();

                            double phi_LH1 = ddvcs_kine.GetPhi_LH_1();
                            double phi_LH2 = ddvcs_kine.GetPhi_LH_2();

                            h_xi_xxGPD_1_4.Fill( xx_GPD1, xi1);
                            h_xi_xxGPD_2_4.Fill( xx_GPD2, xi2);

                            h_Phi_LH1_4.Fill( phi_LH1);
                            h_Phi_LH2_4.Fill( phi_LH2);

                            h_Minv12_4.Fill( ddvcs_kine.GetMinv_1(), ddvcs_kine.GetMinv_2());

                            double dP_P = (Pmis - part_p.p())/part_p.p();
                            h_Pmiss_Pprot4.Fill(part_p.p(), Pmis);
                            h_dP_P_prot4.Fill(part_p.p(), dP_P);

                            if ( TMath::Abs(dP_P) < dP_P_cut ) {
                                h_Mmis6.Fill(Mmis);
                            }

                            if ( TMath::Abs(ddvcs_kine.GetMx2_Reaction() ) < Mx2ReactionCut ) {

                                h_Mmis5.Fill(Mmis);

                                h_xi_xxGPD_1_5.Fill( xx_GPD1, xi1);
                                h_xi_xxGPD_2_5.Fill( xx_GPD2, xi2);
                                h_Minv12_5.Fill( ddvcs_kine.GetMinv_1(), ddvcs_kine.GetMinv_2());
                                h_Phi_LH1_5.Fill( phi_LH1);
                                h_Phi_LH2_5.Fill( phi_LH2);

                                if ( ddvcs_kine.GetMinv_1() > MinvMinCut ) {
                                    h_Phi_LH1_6.Fill( phi_LH1);
                                    h_Phi_LH1_helWeighted_6.Fill( phi_LH1, helicity/beamPol);
                                    h_xi_xxGPD_1_6.Fill( xx_GPD1, xi1);
                                }
                                if ( ddvcs_kine.GetMinv_2() > MinvMinCut ) {
                                    h_Phi_LH2_6.Fill( phi_LH2);
                                    h_Phi_LH2_helWeighted_6.Fill( phi_LH2, helicity/beamPol);
                                    h_xi_xxGPD_2_6.Fill( xx_GPD2, xi2);
                                }
                            }

                        }
                    }

                    if (n_prot == 1) {
                        h_Mmis2.Fill(L_mis.M());
                        h_Mmis_PMis2.Fill(Pmis, Mmis);
                        h_th_VS_Mmis2.Fill(th_mis, Mmis);
                        h_th_P_Mis2.Fill(Pmis, th_mis);

                    }
                }
            }


            if (isMC) {
                event.getStructure(bMCPart);

                double px_ep, py_ep, pz_ep, pt_ep, p_ep;
                double px_em[2], py_em[2], pz_em[2], pt_em[2], p_em[2];

                int pid_em1 = bMCPart.getInt("pid", 1);
                int pid_em2 = bMCPart.getInt("pid", 3);
                int pid_ep = bMCPart.getInt("pid", 2);

                if (pid_em1 != 11 || pid_em2 != 11 || pid_ep != -11) {
                    cout << "This should not happen. Change the code structure!!!!" << endl;
                }

                px_ep = bMCPart.getFloat("px", 2);
                py_ep = bMCPart.getFloat("py", 2);
                pz_ep = bMCPart.getFloat("pz", 2);
                pt_ep = sqrt(px_ep * px_ep + py_ep * py_ep);
                p_ep = sqrt(px_ep * px_ep + py_ep * py_ep + pz_ep * pz_ep);

                px_em[0] = bMCPart.getFloat("px", 1);
                py_em[0] = bMCPart.getFloat("py", 1);
                pz_em[0] = bMCPart.getFloat("pz", 1);
                pt_em[0] = sqrt(px_em[0] * px_em[0] + py_em[0] * py_em[0]);
                p_em[0] = sqrt(px_em[0] * px_em[0] + py_em[0] * py_em[0] + pz_em[0] * pz_em[0]);

                px_em[1] = bMCPart.getFloat("px", 3);
                py_em[1] = bMCPart.getFloat("py", 3);
                pz_em[1] = bMCPart.getFloat("pz", 3);
                pt_em[1] = sqrt(px_em[1] * px_em[1] + py_em[1] * py_em[1]);
                p_em[1] = sqrt(px_em[1] * px_em[1] + py_em[1] * py_em[1] + pz_em[1] * pz_em[1]);

                TLorentzVector L_MC_em1, L_MC_em2, L_MC_ep;
                L_MC_ep.SetPxPyPzE(px_ep, py_ep, pz_ep, p_ep);
                L_MC_em1.SetPxPyPzE(px_em[0], py_em[0], pz_em[0], p_em[0]);
                L_MC_em2.SetPxPyPzE(px_em[1], py_em[1], pz_em[1], p_em[1]);

                bool epAcc = mupAcc(L_MC_ep);
                bool em1_emAcc = emAcc(L_MC_em1);
                bool em1_mum = mumAcc(L_MC_em1);
                bool em2_emAcc = emAcc(L_MC_em2);
                bool em2_mum = mumAcc(L_MC_em2);

                if (epAcc && ((em1_emAcc && em2_mum) || (em2_emAcc && em1_mum))) {
                    TLorentzVector L_emep1 = L_MC_em1 + L_MC_ep;
                    TLorentzVector L_emep2 = L_MC_em2 + L_MC_ep;

                    h_MC_Memep_12_1.Fill(L_emep1.M(), L_emep2.M());
                }
            }
        }
    } catch (const char *msg) {
        cerr << msg << endl;
    }

    /*
        Depending on whether data or MC the scale will be calculated differently,
        In any case it should bring the distribution to "Accepted" cross-section
    */

    double scale = 0.;

    if (isMC) {
        cout << "evCounter = " << evCounter << endl;
        scale = x_sec/evCounter;
    }else {
        scale = 1/Lumi;
    }

    TList *l_objs = gDirectory->GetList();

    for (TIter curObj = l_objs->begin(); curObj != l_objs->end(); curObj.Next()) {
        ((TH1*) * curObj)->Scale(scale);
    }

    /*
     * In order to get into the asymmetry we should divide "hel/pol" weighted histo to non weighted onec
     */
    h_Phi_LH1_helWeighted_6.Divide( &h_Phi_LH1_6 );
    h_Phi_LH2_helWeighted_6.Divide( &h_Phi_LH2_6 );
    h_Phi_LH1_helWeighted_NoProt3.Divide( &h_Phi_LH1_NoProt3 );
    h_Phi_LH2_helWeighted_NoProt3.Divide( &h_Phi_LH2_NoProt3 );


    gDirectory->Write();
    return 0;
}

bool emAcc(TLorentzVector &L) {
    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() <
           P_max;
}

bool mumAcc(TLorentzVector &L) {
    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.5; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() <
           P_max &&
           L.Theta() * TMath::RadToDeg() > f_mumAccThP->Eval(L.P());
}

bool mupAcc(TLorentzVector &L) {
    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 11.; // GeV
    const double P_min = 1.5; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() <
           P_max;
}

bool protAcc(TLorentzVector &L) {
    const double th_max = 120; // deg
    const double th_min = 40.; // deg
    const double P_max = 11.; // GeV
    const double P_min = 0.3; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() <
           P_max;
}
