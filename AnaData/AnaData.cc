/* 
 * File:   AnaData.cc
 * Author: rafopar
 *
 * Created on August 10, 2023, 12:22 PM
 */

#include <cstdlib>
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

#include <RecParticle.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
#include <vector>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {


    char inputFile[256];
    char outputFile[256];
    int run = -9999;
    int fnum = -9999;

    if (argc > 1) {
        run = atoi(argv[1]);
        sprintf(inputFile, "inpHipoFiles/skim4_0%d.hipo", run);
        //sprintf(inputFile, "Data_F18InEarly/jpsitcs_00%d.hipo", run);
        //sprintf(inputFile, "%s", argv[1]);
        sprintf(outputFile, "Skims_RGK_Spring2024/Skim_ememep_%d.hipo", run);
    } else {
        std::cout << " *** please provide the run number..." << std::endl;
        cout << "Exiting ..." << endl;
        exit(0);
    }

    const double Hist_Pmax = 11.;
    const int nMax_samePID = 20;
    const double Mp = 0.9383;
    const double chi2PID_Prot_Max = 10.;
    const double chi2PID_Prot_Min = -10.;
    const int HTCC_TYPE = 15;
    const double HTCC_MomThr = 4.8;
    const double PCalEmin = 0.07;
    const double PCalVmin = 0.07;
    const double PCalWmin = 0.07;
    const double vzMax = 6.;
    const double vzMin = -14.;

    TFile *file_out = new TFile(Form("Skims_RGK_Spring2024/Hists_DDVCS_%d.root", run), "Recreate");
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

    TH2D h_beta_P_prot1("h_beta_P_prot1", "", 200, 0., 1.2 * Hist_Pmax, 200, 0., 1.5);
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



    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;

    hipo::writer writer;
    writer.addDictionary(factory);
    writer.open(outputFile);

    int evCounter = 0;

    hipo::bank bRunConf(factory.getSchema("RUN::config"));
    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));
    hipo::bank bRecSC(factory.getSchema("REC::Scintillator"));
    hipo::bank bRecEV(factory.getSchema("REC::Event"));


    int ind_em[nMax_samePID];
    int ind_ep[nMax_samePID];
    int ind_prot[nMax_samePID];
    int sec_em[nMax_samePID];
    int sec_ep[nMax_samePID];
    double beamCharge = 0.;

    try {

        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            //if( evCounter > 155000 ){break;}
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

            int nPart = bRecPart.getRows();
            for (int i_part = 0; i_part < nPart; i_part++) {

                // ==== Before assigning, index, all indexes are initialized to -1, this way we can check, whether
                // ==== that particular detector is related to the particle "i_part"
                ind_HTCC[i_part] = -1;
                ind_PCal[i_part] = -1;
                ind_ECin[i_part] = -1;
                ind_ECout[i_part] = -1;

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
                RecParticle recp(bRecPart, bRecCalo, bRecCC, ipart, ind_PCal[ipart], ind_ECin[ipart], ind_ECout[ipart], ind_HTCC[ipart]);

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

                    if (((recp.energyECin() < 0.001 && recp.SFPCal() > 0.11) || (recp.energyECin() >= 0.001 && (recp.SFPCal() + recp.SFECin()) > 0.17)) && isPCalEmin && isPCalVmin && isPCalWmin) {
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

            h_n_emep1.Fill(n_em, n_ep);

            if (n_em == 2 && n_ep == 1) {

                h_n_prot1.Fill(n_prot);

                writer.addEvent(event);
                //cout << "Event is is added" << endl;
            }

        }
    } catch (const char msg) {
        cerr << msg << endl;
    }


    gDirectory->Write();
    writer.close();
    writer.showSummary();

    file_out->Close();

    return 0;
}
