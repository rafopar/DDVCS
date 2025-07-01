/* 
 * File:   AnaPionRejection.cc
 * Author: rafopar
 *
 * Created on January 29, 2025, 9:06 AM
 */

#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <RecParticle.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>

using namespace std;

/*
 * This analysis code analyse data with single pion MC events distributed uniformly in momentum and theta and phi,
 * then will study what fraction of pions are detected as muons.
 */


struct MCPart {
    int pid;
    double px;
    double py;
    double pz;
    double p;
    double vx;
    double vy;
    double vz;
    double theta;
    double phi;
};

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Please provide the Identifier as an argument" << endl;
        cout << "The program should run as ./AnaPionRejection.exe Identifier" << endl;
        cout << "As an example ./AnaPionRejection.exe piPlus_6GeV_22Deg_0Deg_5GeV_17deg_40deg" << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    std::string Identifier = argv[1];

    cout << "Identifier is " << Identifier << endl;

    const int HTCC_TYPE = 15;
    const double Emax = 12.;
    const double cosThMax = cos(41 * TMath::DegToRad());
    const double d2r = TMath::DegToRad();

    const int nMaxPCalStripsPos = 5;
    const int nMaxECinStripsPos = 5;
    const int nMaxECoutStripsPos = 5;
    const int nMaxPCalStripsNeg = 7;
    const int nMaxECinStripsNeg = 5;
    const int nMaxECoutStripsNeg = 5;

    const double LoosPCal_Ecut_Pos = 0.075;
    const double LoosECin_Ecut_Pos = 0.075;
    const double LoosECout_Ecut_Pos = 0.094;

    const double LoosECout_Ecut_Neg = 0.1;
    const double LoosECin_Ecut_Neg = 0.075;
    const double LoosPCal_Ecut_Neg = 0.068;

    TFile *file_out = new TFile(Form("pion_Rejection_%s.root", Identifier.c_str()), "Recreate");

    TH1D h_pos_E_PCal1("h_pos_E_PCal1", "", 200, 0., 0.25);
    TH1D h_pos_E_PCal2("h_pos_E_PCal2", "", 200, 0., 0.25);
    TH1D h_pos_E_ECin1("h_pos_E_ECin1", "", 200, 0., 0.25);
    TH1D h_pos_E_ECin2("h_pos_E_ECin2", "", 200, 0., 0.25);
    TH1D h_pos_E_ECout1("h_pos_E_ECout1", "", 200, 0., 0.25);
    TH1D h_pos_E_ECout2("h_pos_E_ECout2", "", 200, 0., 0.25);

    TH1D h_neg_E_PCal1("h_neg_E_PCal1", "", 200, 0., 0.25);
    TH1D h_neg_E_PCal2("h_neg_E_PCal2", "", 200, 0., 0.25);
    TH1D h_neg_E_ECin1("h_neg_E_ECin1", "", 200, 0., 0.25);
    TH1D h_neg_E_ECin2("h_neg_E_ECin2", "", 200, 0., 0.25);
    TH1D h_neg_E_ECout1("h_neg_E_ECout1", "", 200, 0., 0.25);
    TH1D h_neg_E_ECout2("h_neg_E_ECout2", "", 200, 0., 0.25);

    TH2D h_SF_pos1("h_SF_pos1", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_pos2("h_SF_pos2", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_pos3("h_SF_pos3", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_pos4("h_SF_pos4", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_pos5("h_SF_pos5", "", 200, 0., Emax, 200, 0., 0.35);

    TH2D h_SF_neg1("h_SF_neg1", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_neg2("h_SF_neg2", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_neg3("h_SF_neg3", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_neg4("h_SF_neg4", "", 200, 0., Emax, 200, 0., 0.35);
    TH2D h_SF_neg5("h_SF_neg5", "", 200, 0., Emax, 200, 0., 0.35);

    TH1D h_pos_n_StripPCal1("h_pos_n_StripPCal1", "", 11, -0.5, 10.5);
    TH1D h_pos_n_StripPCal2("h_pos_n_StripPCal2", "", 11, -0.5, 10.5);
    TH1D h_pos_n_StripECin1("h_pos_n_StripECin1", "", 11, -0.5, 10.5);
    TH1D h_pos_n_StripECin2("h_pos_n_StripECin2", "", 11, -0.5, 10.5);
    TH1D h_pos_n_StripECout1("h_pos_n_StripECout1", "", 11, -0.5, 10.5);
    TH1D h_pos_n_StripECout2("h_pos_n_StripECout2", "", 11, -0.5, 10.5);

    TH1D h_neg_n_StripPCal1("h_neg_n_StripPCal1", "", 11, -0.5, 10.5);
    TH1D h_neg_n_StripPCal2("h_neg_n_StripPCal2", "", 11, -0.5, 10.5);
    TH1D h_neg_n_StripECin1("h_neg_n_StripECin1", "", 11, -0.5, 10.5);
    TH1D h_neg_n_StripECin2("h_neg_n_StripECin2", "", 11, -0.5, 10.5);
    TH1D h_neg_n_StripECout1("h_neg_n_StripECout1", "", 11, -0.5, 10.5);
    TH1D h_neg_n_StripECout2("h_neg_n_StripECout2", "", 11, -0.5, 10.5);

    TH2D h_cosTh_P_pi_MC1("h_cosTh_P_pi_MC1", "", 15, 0., Emax, 15., cosThMax, 1.);
    TH2D h_cosTh_P_pi_mup_MC2("h_cosTh_P_pi_mup_MC2", "", 15, 0., Emax, 15., cosThMax, 1.);
    TH2D h_cosTh_P_pi_mup_MC3("h_cosTh_P_pi_mup_MC3", "", 15, 0., Emax, 15., cosThMax, 1.);
    TH2D h_cosTh_P_pi_mum_MC2("h_cosTh_P_pi_mum_MC2", "", 15, 0., Emax, 15., cosThMax, 1.);
    TH2D h_cosTh_P_pi_mum_MC3("h_cosTh_P_pi_mum_MC3", "", 15, 0., Emax, 15., cosThMax, 1.);
    TH2D h_cosTh_phi_pi_MC1("h_cosTh_phi_pi_MC1", "", 50, -45., 45., 50., cosThMax, 1.);
    TH2D h_cosTh_phi_pi_mup_MC2("h_cosTh_phi_pi_mup_MC2", "", 50, -45., 45., 50., cosThMax, 1.);
    TH2D h_cosTh_phi_pi_mup_MC3("h_cosTh_phi_pi_mup_MC3", "", 50, -45., 45., 50., cosThMax, 1.);
    TH2D h_cosTh_phi_pi_mum_MC2("h_cosTh_phi_pi_mum_MC2", "", 50, -45., 45., 50., cosThMax, 1.);
    TH2D h_cosTh_phi_pi_mum_MC3("h_cosTh_phi_pi_mum_MC3", "", 50, -45., 45., 50., cosThMax, 1.);

    TH1D h_vz_pi_MC1("h_vz_pi_MC1", "", 1200, -6., 6.);
    TH1D h_vy_pi_MC1("h_vy_pi_MC1", "", 1200, -6., 6.);
    TH1D h_vx_pi_MC1("h_vx_pi_MC1", "", 1200, -6., 6.);
    
    TH2D h_cosTh_P_pi_mup_Loose_MC2("h_cosTh_P_pi_mup_Loose_MC2", "", 15, 0., Emax, 15., cosThMax, 1.);
    TH2D h_cosTh_P_pi_mum_Loose_MC2("h_cosTh_P_pi_mum_Loose_MC2", "", 15, 0., Emax, 15., cosThMax, 1.);

    TH2D h_dP_P_mup1("h_dP_P_mup1", "", 50, 0., Emax, 50, -1., 0.05);
    TH2D h_dP_P_mum1("h_dP_P_mum1", "", 50, 0., Emax, 50, -1., 0.05);

    TH2D h_PRec_PMC_mup1("h_PRec_PMC_mup1", "", 50, 0., Emax, 50, 0., Emax);
    TH2D h_PRec_PMC_mum1("h_PRec_PMC_mum1", "", 50, 0., Emax, 50, 0., Emax);

    TH1D h_N_mup1("h_N_mup1", "", 5, -0.5, 4.5);
    TH1D h_N_mum1("h_N_mum1", "", 5, -0.5, 4.5);

    TH1D h_P_mup1("h_P_mup1", "", 200, 0., Emax);
    TH1D h_P_mup2("h_P_mup2", "", 200, 0., Emax);
    TH1D h_P_mum1("h_P_mum1", "", 200, 0., Emax);
    TH1D h_P_mum2("h_P_mum2", "", 200, 0., Emax);

    hipo::reader reader;
    //reader.open("Data/Skim_ZeroSuppr_2247_All.hipo");
    reader.open(Form("Data/%s_All.hipo", Identifier.c_str()));

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();
    hipo::event event;
    int evCounter = 0;

    hipo::bank bMCPart(factory.getSchema("MC::Particle"));
    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));
    hipo::bank bRecSC(factory.getSchema("REC::Scintillator"));
    hipo::bank bRecEV(factory.getSchema("REC::Event"));
    hipo::bank bRunConf(factory.getSchema("RUN::config"));

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

            event.getStructure(bMCPart);
            event.getStructure(bRecPart);
            event.getStructure(bRecCalo);

            map<int, int> ind_HTCC;
            map<int, int> ind_PCal;
            map<int, int> ind_ECin;
            map<int, int> ind_ECout;

            int nPart = bRecPart.getRows();
            int nMCPart = bMCPart.getRows();


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

            vector<RecParticle> v_recp_mum;
            vector<RecParticle> v_recp_mup;

            vector<RecParticle> v_recp_mum_looseTrg;
            vector<RecParticle> v_recp_mup_looseTrg;

            for (int ipart = 0; ipart < nPart; ipart++) {
                RecParticle recp(bRecPart, bRecCalo, bRecCC, ipart, ind_PCal[ipart], ind_ECin[ipart], ind_ECout[ipart], ind_HTCC[ipart]);


                if ((recp.pid() == 211 || recp.pid() == -13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {
                    h_pos_E_PCal1.Fill(recp.energyPCal());
                    h_pos_E_ECin1.Fill(recp.energyECin());
                    h_pos_E_ECout1.Fill(recp.energyECout());

                    h_SF_pos1.Fill(recp.p(), recp.SF());

                    int nStrip_PCal = recp.duPCal() + recp.dvPCal() + recp.dwPCal();
                    int nStrip_ECin = recp.duECin() + recp.dvECin() + recp.dwECin();
                    int nStrip_ECout = recp.duECout() + recp.dvECout() + recp.dwECout();

                    bool hitAllLayers = (nStrip_PCal >= 3) && (nStrip_ECin >= 3) && (nStrip_ECout >= 3);

                    h_pos_n_StripPCal1.Fill(nStrip_PCal);
                    h_pos_n_StripECin1.Fill(nStrip_ECin);
                    h_pos_n_StripECout1.Fill(nStrip_ECout);

                    if (recp.energyPCal() > 0) {
                        h_SF_pos2.Fill(recp.p(), recp.SF());
                    }
                    if (recp.energyECin() > 0) {
                        h_SF_pos3.Fill(recp.p(), recp.SF());
                    }
                    if (recp.energyECout() > 0 && recp.energyPCal() > 0 && recp.energyECin() > 0) {
                        h_SF_pos4.Fill(recp.p(), recp.SF());
                        h_pos_E_PCal2.Fill(recp.energyPCal());
                        h_pos_E_ECin2.Fill(recp.energyECin());
                        h_pos_E_ECout2.Fill(recp.energyECout());

                        h_pos_n_StripPCal2.Fill(nStrip_PCal);
                        h_pos_n_StripECin2.Fill(nStrip_ECin);
                        h_pos_n_StripECout2.Fill(nStrip_ECout);

                        /* 
                         * Identifying muons with a loose energy cuts that will be applied on trigger
                         */

                        if (recp.energyECout() < LoosECout_Ecut_Pos && recp.energyECin() < LoosECin_Ecut_Pos && recp.energyPCal() < LoosPCal_Ecut_Pos) {
                            v_recp_mup_looseTrg.push_back(recp);
                        }

                        /* 
                         * Identifying muons with a tighter cuts, like in the analysis ti suppress pions
                         */

                        if (nStrip_PCal <= nMaxPCalStripsPos && nStrip_ECin < nMaxECinStripsPos && nStrip_ECout < nMaxECoutStripsPos && hitAllLayers) {
                            h_SF_pos5.Fill(recp.p(), recp.SF());
                            v_recp_mup.push_back(recp);
                        }
                    }
                } else if ((recp.pid() == -211 || recp.pid() == 13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {

                    h_neg_E_PCal1.Fill(recp.energyPCal());
                    h_neg_E_ECin1.Fill(recp.energyECin());
                    h_neg_E_ECout1.Fill(recp.energyECout());

                    h_SF_neg1.Fill(recp.p(), recp.SF());

                    int nStrip_PCal = recp.duPCal() + recp.dvPCal() + recp.dwPCal();
                    int nStrip_ECin = recp.duECin() + recp.dvECin() + recp.dwECin();
                    int nStrip_ECout = recp.duECout() + recp.dvECout() + recp.dwECout();
                    bool hitAllLayers = (nStrip_PCal >= 3) && (nStrip_ECin >= 3) && (nStrip_ECout >= 3);

                    h_neg_n_StripPCal1.Fill(nStrip_PCal);
                    h_neg_n_StripECin1.Fill(nStrip_ECin);
                    h_neg_n_StripECout1.Fill(nStrip_ECout);

                    if (recp.energyPCal() > 0) {
                        h_SF_neg2.Fill(recp.p(), recp.SF());
                    }
                    if (recp.energyECin() > 0) {
                        h_SF_neg3.Fill(recp.p(), recp.SF());
                    }
                    if (recp.energyECout() > 0 && recp.energyPCal() > 0 && recp.energyECin() > 0) {
                        h_SF_neg4.Fill(recp.p(), recp.SF());
                        h_neg_E_PCal2.Fill(recp.energyPCal());
                        h_neg_E_ECin2.Fill(recp.energyECin());
                        h_neg_E_ECout2.Fill(recp.energyECout());

                        h_neg_n_StripPCal2.Fill(nStrip_PCal);
                        h_neg_n_StripECin2.Fill(nStrip_ECin);
                        h_neg_n_StripECout2.Fill(nStrip_ECout);

                        /* 
                         * Identifying muons with a loose energy cuts that will be applied on trigger
                         */

                        if (recp.energyECout() < LoosECout_Ecut_Neg && recp.energyECin() < LoosECin_Ecut_Neg && recp.energyPCal() < LoosPCal_Ecut_Neg) {
                            v_recp_mum_looseTrg.push_back(recp);
                        }

                        /* 
                         * Identifying muons with a tighter cuts, like in the analysis ti suppress pions
                         */



                        if (nStrip_PCal <= nMaxPCalStripsNeg && nStrip_ECin < nMaxECinStripsNeg && nStrip_ECout < nMaxECoutStripsNeg && hitAllLayers) {
                            h_SF_neg5.Fill(recp.p(), recp.SF());

                            v_recp_mum.push_back(recp);
                        }
                    }
                }
            }

            MCPart mc_pi;

            mc_pi.pid = bMCPart.getInt("pid", 0);
            mc_pi.px = double(bMCPart.getFloat("px", 0));
            mc_pi.py = double(bMCPart.getFloat("py", 0));
            mc_pi.pz = double(bMCPart.getFloat("pz", 0));
            mc_pi.p = sqrt(mc_pi.px * mc_pi.px + mc_pi.py * mc_pi.py + mc_pi.pz * mc_pi.pz);
            mc_pi.vx = double(bMCPart.getFloat("vx", 0));
            mc_pi.vy = double(bMCPart.getFloat("vy", 0));
            mc_pi.vz = double(bMCPart.getFloat("vz", 0));
            mc_pi.theta = atan(sqrt(mc_pi.px * mc_pi.px + mc_pi.py * mc_pi.py) / mc_pi.pz) * TMath::RadToDeg();
            mc_pi.phi = atan2(mc_pi.py, mc_pi.px) * TMath::RadToDeg();

            h_vz_pi_MC1.Fill(mc_pi.vz);
            h_vy_pi_MC1.Fill(mc_pi.vy);
            h_vx_pi_MC1.Fill(mc_pi.vx);
            
            h_cosTh_P_pi_MC1.Fill(mc_pi.p, cos(mc_pi.theta * d2r));
            h_cosTh_phi_pi_MC1.Fill(mc_pi.phi, cos(mc_pi.theta * d2r));

            int n_mup = v_recp_mup.size();
            int n_mum = v_recp_mum.size();

            h_N_mup1.Fill(n_mup);
            h_N_mum1.Fill(n_mum);

            for (auto curP : v_recp_mup) {
                h_P_mup1.Fill(curP.p());
                h_dP_P_mup1.Fill(curP.p(), (curP.p() - mc_pi.p) / mc_pi.p);
                h_PRec_PMC_mup1.Fill(mc_pi.p, curP.p());

                if (curP.p() > -2.5 + mc_pi.p && curP.p() < -0.8 + mc_pi.p) {
                    h_P_mup2.Fill(curP.p());
                    h_cosTh_P_pi_mup_MC3.Fill(mc_pi.p, cos(mc_pi.theta * d2r));
                }

            }

            for (auto curP : v_recp_mum) {
                h_P_mum1.Fill(curP.p());
                h_dP_P_mum1.Fill(curP.p(), (curP.p() - mc_pi.p) / mc_pi.p);
                h_PRec_PMC_mum1.Fill(mc_pi.p, curP.p());

                if (curP.p() > -2.5 + mc_pi.p && curP.p() < -0.8 + mc_pi.p) {
                    h_P_mum2.Fill(curP.p());
                }
            }

            if (v_recp_mup_looseTrg.size() > 0) {
                h_cosTh_P_pi_mup_Loose_MC2.Fill(mc_pi.p, cos(mc_pi.theta * d2r));
            }

            if (v_recp_mum_looseTrg.size() > 0) {
                h_cosTh_P_pi_mum_Loose_MC2.Fill(mc_pi.p, cos(mc_pi.theta * d2r));
            }

            if (v_recp_mup.size() > 0) {
                h_cosTh_P_pi_mup_MC2.Fill(mc_pi.p, cos(mc_pi.theta * d2r));
                h_cosTh_phi_pi_mup_MC2.Fill(mc_pi.phi, cos(mc_pi.theta * d2r));
            }

            if (v_recp_mum.size() > 0) {
                h_cosTh_P_pi_mum_MC2.Fill(mc_pi.p, cos(mc_pi.theta * d2r));
                h_cosTh_phi_pi_mum_MC2.Fill(mc_pi.phi, cos(mc_pi.theta * d2r));
            }

        }
    } catch (const char *msg) {
        cerr << msg << endl;
    }

    gDirectory->Write();

    return 0;
}