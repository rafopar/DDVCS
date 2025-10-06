/* 
 * File:   SkipPionPair.cc
 * Author: rafopar
 *
 * Created on March 7, 2025, 10:42 AM
 */

#include <cstdlib>
#include <fstream>

#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TRandom.h>
#include <TLorentzVector.h>

#include <RecParticle.h>

// ===== Hipo headers =====
#include <reader.h>
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 4) {
        cout << "Please provide the run number and file Identifier N_PerLUND" << endl;
        cout << "The command should look like ./SkipPionPair.exe 5117 00000-00004 1000" << endl;
        cout << "In the above example, the run is 5117m file Identifier is \"00000-00004\" and N_PerLUND = 1000 " << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    int run = atoi(argv[1]);
    std::string file_ID = argv[2];
    int NeVPerFile = atoi(argv[3]);

    const int HTCC_TYPE = 15;
    const double Mmu = 0.10566;
    const double Mpi = 0.13957039;
    const double Mp = 0.9383;
    const double Eb = 10.6;
    const double r2d = TMath::RadToDeg();
    const double th_MuMax = 40; // We ignore angles above this value
    const bool writeLUND = true;
    const bool dumpEvNumbers = false;

    const double PMinCut = 2;
    const double Mmis_max = 1.2;
    const double Mmis_min = 0.6;
    const double LowMinvCut = 1.; // We want to study those evens with Minv < 1 GeV
    const int nMaxPCalStripsPos = 5;
    const int nMaxECinStripsPos = 5;
    const int nMaxECoutStripsPos = 5;
    const int nMaxPCalStripsNeg = 7;
    const int nMaxECinStripsNeg = 5;
    const int nMaxECoutStripsNeg = 5;

    TLorentzVector L_beam(0., 0., Eb, Eb);
    TLorentzVector L_targ(0., 0., 0., Mp);

    TFile *file_out = new TFile(Form("Skim_PionPair_%d_%s.root", run, file_ID.c_str()), "Recreate");

    TH1D h_nPip_MIP1("h_nPip_MIP1", "", 5, -0.5, 4.5);
    TH1D h_nPim_MIP1("h_nPim_MIP1", "", 5, -0.5, 4.5);
    TH1D h_Minv_pion_MIP1("h_Minv_pion_MIP1", "", 200, 0., 4);
    TH1D h_Minv_pion_MIP2("h_Minv_pion_MIP2", "", 200, 0., 4);
    TH1D h_Minv_pion_MIP3("h_Minv_pion_MIP3", "", 200, 0., 4);
    TH1D h_Minv_pion_MIP4("h_Minv_pion_MIP4", "", 200, 0., 4);
    TH1D h_Minv_pion_MIP5("h_Minv_pion_MIP5", "", 200, 0., 4);
    TH1D h_Minv_pion_MIP6("h_Minv_pion_MIP6", "", 200, 0., 4);
    TH1D h_Mmis1("h_Mmis1", "", 200, 0., 3.);
    TH1D h_Mmis2("h_Mmis2", "", 200, 0., 3.);
    TH2D h_th_P_pip1("h_th_P_pip1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_pip2("h_th_P_pip2", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_pim1("h_th_P_pim1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_pim2("h_th_P_pim2", "", 200, 0., 1.05 * Eb, 200, 0., 45);

    hipo::reader reader;
    //reader.open("Data/Skim_ZeroSuppr_2247_All.hipo");
    reader.open(Form("Data/%d/rec_clas_00%d_%s.hipo", run, run, file_ID.c_str()));

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();
    hipo::event event;
    int evCounter = 0;

    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));
    hipo::bank bRecSC(factory.getSchema("REC::Scintillator"));
    hipo::bank bRecEV(factory.getSchema("REC::Event"));
    hipo::bank bRunConf(factory.getSchema("RUN::config"));

    int LUND_Fileindex = 0;
    int nLUND_Event = 0;
    int nTot_Event = 0;
    ofstream lund_Out;
    ofstream out_EvNumbers;

    if (dumpEvNumbers) {
        out_EvNumbers.open(Form("EvNumbers_Run_%d_files_%s.dat", run, file_ID.c_str()));
    }

    if (writeLUND) {
        lund_Out.open(Form("Out_LUND/%d/pion_PairSkim_%d_hipoInd_%s_LUNDIndex_%d.txt", run, run, file_ID.c_str(), LUND_Fileindex));
    }

    try {
        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            //            if (evCounter > 200000) {
            //                break;
            //            }
            if (evCounter % 10000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

            if (writeLUND) {
                //if ((nLUND_Event + 1) % NeVPerFile == 0) {
                if (nLUND_Event >= NeVPerFile) {
                    lund_Out.close();
                    LUND_Fileindex++;
                    lund_Out.open(Form("Out_LUND/%d/pion_PairSkim_%d_hipoInd_%s_LUNDIndex_%d.txt", run, run, file_ID.c_str(), LUND_Fileindex));
                    nLUND_Event = 0;
                    nTot_Event++;
                }
            }

            event.getStructure(bRecPart);
            event.getStructure(bRecCalo);
            event.getStructure(bRunConf);

            int ev_Number = bRunConf.getInt("event", 0);

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

            vector<RecParticle> v_recp_pim, v_recp_pim_MIP;
            vector<RecParticle> v_recp_pip, v_recp_pip_MIP;
            vector<RecParticle> v_recp_em;

            for (int ipart = 0; ipart < nPart; ipart++) {
                RecParticle recp(bRecPart, bRecCalo, bRecCC, ipart, ind_PCal[ipart], ind_ECin[ipart], ind_ECout[ipart], ind_HTCC[ipart]);

                if ((recp.pid() == 211 || recp.pid() == -13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {
                    v_recp_pip.push_back(recp);

                    int nStrip_PCal = recp.duPCal() + recp.dvPCal() + recp.dwPCal();
                    int nStrip_ECin = recp.duECin() + recp.dvECin() + recp.dwECin();
                    int nStrip_ECout = recp.duECout() + recp.dvECout() + recp.dwECout();

                    bool hitAllLayers = (nStrip_PCal >= 3) && (nStrip_ECin >= 3) && (nStrip_ECout >= 3);

                    if (recp.th() < th_MuMax) {
                        bool MIP_Signature = nStrip_PCal <= nMaxPCalStripsPos && nStrip_ECin < nMaxECinStripsPos && nStrip_ECout < nMaxECoutStripsPos;

                        if (!MIP_Signature && hitAllLayers) {
                            v_recp_pip_MIP.push_back(recp);
                            continue;
                        }
                    }

                } else if ((recp.pid() == -211 || recp.pid() == 13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {
                    v_recp_pim.push_back(recp);

                    int nStrip_PCal = recp.duPCal() + recp.dvPCal() + recp.dwPCal();
                    int nStrip_ECin = recp.duECin() + recp.dvECin() + recp.dwECin();
                    int nStrip_ECout = recp.duECout() + recp.dvECout() + recp.dwECout();

                    bool hitAllLayers = (nStrip_PCal >= 3) && (nStrip_ECin >= 3) && (nStrip_ECout >= 3);

                    if (recp.th() < th_MuMax) {
                        bool MIP_Signature = nStrip_PCal <= nMaxPCalStripsNeg && nStrip_ECin < nMaxECinStripsNeg && nStrip_ECout < nMaxECoutStripsNeg;

                        if (!MIP_Signature && hitAllLayers) {
                            v_recp_pim_MIP.push_back(recp);
                            continue;
                        }
                    }

                } else if (recp.pid() == 11) {
                    v_recp_em.push_back(recp);
                }
            }


            h_nPip_MIP1.Fill(v_recp_pip_MIP.size());
            h_nPim_MIP1.Fill(v_recp_pim_MIP.size());

            if (v_recp_pim_MIP.size() >= 1 && v_recp_pip_MIP.size() >= 1) {

                TLorentzVector L_pip_MIP, L_pim_MIP, L_mis, L_em;

                for (auto rec_pip : v_recp_pip_MIP) {

                    h_th_P_pip1.Fill(rec_pip.p(), rec_pip.th());

                    for (auto rec_pim : v_recp_pim_MIP) {

                        L_pip_MIP.SetPxPyPzE(rec_pip.px(), rec_pip.py(), rec_pip.pz(), sqrt(rec_pip.p() * rec_pip.p() + Mpi * Mpi));
                        L_pim_MIP.SetPxPyPzE(rec_pim.px(), rec_pim.py(), rec_pim.pz(), sqrt(rec_pim.p() * rec_pim.p() + Mpi * Mpi));

                        TLorentzVector L_pimpip = L_pim_MIP + L_pip_MIP;
                        double Minv = L_pimpip.M();
                        h_Minv_pion_MIP1.Fill(Minv);

                        double Mx2;

                        //Mmis_max
                        bool MissMassCut;

                        if (v_recp_em.size() >= 1) {
                            h_Minv_pion_MIP3.Fill(Minv);

                            L_em.SetPxPyPzE(v_recp_em.at(0).px(), v_recp_em.at(0).py(), v_recp_em.at(0).pz(), v_recp_em.at(0).p());

                            L_mis = L_beam + L_targ - L_em - L_pimpip;

                            Mx2 = L_mis.M2();

                            MissMassCut = Mx2 > Mmis_min && Mx2 < Mmis_max;
                            h_Mmis1.Fill(Mx2);

                            if (MissMassCut) {
                                h_Minv_pion_MIP5.Fill(Minv);
                            }
                        }

                        if (L_pip_MIP.P() > PMinCut && L_pim_MIP.P() > PMinCut) {
                            h_Minv_pion_MIP2.Fill(Minv);
                            if (v_recp_em.size() >= 1) {
                                h_Minv_pion_MIP4.Fill(Minv);
                                h_Mmis2.Fill(Mx2);

                                if (MissMassCut) {
                                    h_Minv_pion_MIP6.Fill(Minv);

                                    if (Minv < LowMinvCut) {
                                        h_th_P_pip2.Fill(rec_pip.p(), rec_pip.th());
                                        h_th_P_pim2.Fill(rec_pim.p(), rec_pim.th());
                                    }

                                }

                            }
                        }
                    }
                }

                for (auto rec_pim : v_recp_pim_MIP) {
                    h_th_P_pim1.Fill(rec_pim.p(), rec_pim.th());
                }
            }

            /* 
             * Now let's check if there is at least one (pim,pip) pair, then write this to a LUND file
             */

            if (!(v_recp_pip.size() >= 1 && v_recp_pim.size() >= 1)) {
                continue;
            }

            out_EvNumbers << ev_Number << endl;

            int n_pip = v_recp_pip.size();
            int n_pim = v_recp_pim.size();
            int n_em = v_recp_em.size();

            int n_Tot = n_em + n_pim + n_pip;

            int UD = 0; // Used Designed, Not used by GEMC
            double beam_Pol = 0.8; // Although dosn't matter

            if (writeLUND) {
                    
                    lund_Out << n_Tot << setw(5) << UD << setw(5) << UD << setw(5) << UD << setw(5) << beam_Pol << setw(5) << UD << setw(5) << UD << setw(5) << UD << setw(15) << ev_Number << setw(5) << UD << endl;

                    int lund_ind = 1;

                    /* 
                     * Writing pi+
                     */
                    for (auto cur_pip : v_recp_pip) {
                        int pid = cur_pip.pid();
                        double vx = 0;
                        double vy = 0;
                        double vz = 0;

                        lund_Out << lund_ind << setw(5) << UD << setw(5) << 1 << setw(5) << pid << setw(5) << UD << setw(5) << UD << setw(15) << cur_pip.px() << setw(15) << cur_pip.py() << setw(15) << cur_pip.pz() << UD << setw(5) << UD <<
                                setw(5) << vx << setw(5) << vy << setw(5) << vz << endl;
                        lund_ind = lund_ind + 1;
                    }


                    /* 
                     * Writing pi-
                     */
                    for (auto cur_pim : v_recp_pim) {
                        int pid = cur_pim.pid();
                        double vx = 0;
                        double vy = 0;
                        double vz = 0;

                        lund_Out << lund_ind << setw(5) << UD << setw(5) << 1 << setw(5) << pid << setw(5) << UD << setw(5) << UD << setw(15) << cur_pim.px() << setw(15) << cur_pim.py() << setw(15) << cur_pim.pz() << UD << setw(5) << UD <<
                                setw(5) << vx << setw(5) << vy << setw(5) << vz << endl;
                        lund_ind = lund_ind + 1;
                    }

                    /* 
                     * Writing e-
                     */
                    for (auto cur_em : v_recp_em) {
                        int pid = cur_em.pid();
                        double vx = 0;
                        double vy = 0;
                        double vz = 0;

                        lund_Out << lund_ind << setw(5) << UD << setw(5) << 1 << setw(5) << pid << setw(5) << UD << setw(5) << UD << setw(15) << cur_em.px() << setw(15) << cur_em.py() << setw(15) << cur_em.pz() << UD << setw(5) << UD <<
                                setw(5) << vx << setw(5) << vy << setw(5) << vz << endl;
                        lund_ind = lund_ind + 1;
                    }
                    nLUND_Event = nLUND_Event + 1;
                    nTot_Event = nTot_Event + 1;

            }
        }
    } catch (const char *msg) {
        cerr << msg << endl;
    }

    cout << endl << endl;
    cout << "Total " << nTot_Event << endl;
    cout << "Run = " << run << endl;
    cout << "file_ID = " << file_ID << endl;
    cout << "NeVPerFile = " << NeVPerFile << endl;

    gDirectory->Write();

    return 0;
}