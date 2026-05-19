//
// Created by rafopar on 12/22/25.
//

#include <iostream>
#include <ostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <reader.h>
#include <dictionary.h>

#include <ElectronDDVCSKine.h>
#include <RecParticle.h>

using namespace std;
using LVec = ROOT::Math::PxPyPzEVector;

int main( int argc, char *argv[] ) {

    char inputFile[256];

    if( argc < 2 ) {
        cerr << "Usage: CreateTrainingSample.exe <Run>" << endl;
        cerr << "Exiting" << endl;
        return EXIT_FAILURE;
    }

    const double Eb = 10.6;
    int run = atoi(argv[1]);
    sprintf(inputFile, "Data/Run%d_TrainSample.hipo", run);

    constexpr double r2d = TMath::RadToDeg();
    constexpr double Mp = 0.9383;
    constexpr double MxRecoil_Max = 1.1; // GeV
    constexpr double MxRecoil_Min = 0.88; // GeV

    ElectronDDVCSKine ddvcs_kine_Noprot;
    ddvcs_kine_Noprot.SetEb(Eb);
    ddvcs_kine_Noprot.SetMtarg(Mp);


    TFile outFile(Form("TrainingSample/TrainingSample_Run_%d.root", run), "Recreate");
    TH1D h_vz_em1("h_vz_em1", "", 200, -15., 15.);
    TH1D h_vz_ep1("h_vz_ep1", "", 200, -15., 15.);
    TH1D h_vz_prot("h_vz_prot", "", 200, -15., 15.);

    TH1D h_PCal_lU_em1("h_PCal_lU_em1", "", 200, 0., 450);
    TH1D h_PCal_lV_em1("h_PCal_lV_em1", "", 200, 0., 450);
    TH1D h_PCal_lW_em1("h_PCal_lW_em1", "", 200, 0., 450);
    TH1D h_PCal_lU_ep1("h_PCal_lU_ep1", "", 200, 0., 450);
    TH1D h_PCal_lV_ep1("h_PCal_lV_ep1", "", 200, 0., 450);
    TH1D h_PCal_lW_ep1("h_PCal_lW_ep1", "", 200, 0., 450);

    TH2D h_quality_em2_em1_1("h_quality_em2_em1_1", "", 200, -0.1, 1.1, 200, -0.1, 1.1);
    TH1D h_MatchQuality_em_beam1("h_MatchQuality_em_beam1", "", 200, -0.1, 1.01);
    TH1D h_MatchQuality_em_beam2("h_MatchQuality_em_beam2", "", 200, -0.1, 1.01);
    TH1D h_MatchQuality_em_decay1("h_MatchQuality_em_decay1", "", 200, -0.1, 1.01);
    TH1D h_MatchQuality_em_decay2("h_MatchQuality_em_decay2", "", 200, -0.1, 1.01);

    TH2D h_ind_embeam_decay1("h_ind_embeam_decay1", "", 11, -1.5, 9.5, 11, -1.5, 9.5);

    TH2D h_th_p_Rec_emDecay1("h_th_p_Rec_emDecay1", "", 200, 0., 11., 200, 0, 40.);
    TH2D h_th_p_Rec_emBeam1("h_th_p_Rec_emBeam1", "", 200, 0., 11., 200, 0, 40.);
    TH2D h_th_p_Rec_ep1("h_th_p_Rec_ep1", "", 200, 0., 11., 200, 0, 40.);

    TH2D h_th_p_MC_emBeam1("h_th_p_MC_emBeam1", "", 200, 0., 11., 200, 0, 40.);
    TH2D h_th_p_MC_emDecay1("h_th_p_MC_emDecay1", "", 200, 0., 11., 200, 0, 40.);
    TH2D h_th_p_MC_ep1("h_th_p_MC_ep1", "", 200, 0., 11., 200, 0, 40.);

    TH2D h_th_p_MC_emBeam2("h_th_p_MC_emBeam2", "", 200, 0., 11., 200, 0, 40.);
    TH2D h_th_p_MC_emDecay2("h_th_p_MC_emDecay2", "", 200, 0., 11., 200, 0, 40.);
    TH2D h_th_p_MC_ep2("h_th_p_MC_ep2", "", 200, 0., 11., 200, 0, 40.);

    TH2D h_P_em_Decay_Beam1("h_P_em_Decay_Beam1", "", 200, 0., 11., 200, 0., 11);
    TH2D h_P_emep_Decay1("h_P_emep_Decay1", "", 200, 0., 11., 200, 0., 11);
    TH2D h_P_emep_Beam1("h_P_emep_Beam1", "", 200, 0., 11., 200, 0., 11);

    TH2D h_Q2_Decay_Beam1("h_Q2_Decay_Beam1", "", 200, 0., 2.5, 200, 0., 2.5);
    TH2D h_Q2_Decay_Beam2("h_Q2_Decay_Beam2", "", 200, 0., 2.5, 200, 0., 2.5);

    TH2D h_Minv_Beam_Decay_MC1("h_Minv_Beam_Decay_MC1", "", 200, 0., 3., 200, 0., 3.);
    TH2D h_Minv_Beam_Decay_MC2("h_Minv_Beam_Decay_MC2", "", 200, 0., 3., 200, 0., 3.);

    /*
     * Histograms for assessing the quality of the Truth Matching
     */

    TH2D h_P_beam_em_Rec_MC1("h_P_beam_em_Rec_MC1", "", 200, 0., 11., 200, 0., 11.);
    TH2D h_th_beam_em_Rec_MC1("h_th_beam_em_Rec_MC1", "", 200., 0., 40., 200, 0., 40);
    TH2D h_phi_beam_em_Rec_MC1("h_phi_beam_em_Rec_MC1", "", 200, 0., 360., 200, 0., 360);
    TH2D h_P_decay_em_Rec_MC1("h_P_decay_em_Rec_MC1", "", 200, 0., 11., 200, 0., 11.);
    TH2D h_th_decay_em_Rec_MC1("h_th_decay_em_Rec_MC1", "", 200., 0., 40., 200, 0., 40);
    TH2D h_phi_decay_em_Rec_MC1("h_phi_decay_em_Rec_MC1", "", 200, 0., 360., 200, 0., 360);

    TH1D h_Mmis1("h_Mmis1", "", 200, 0., 4.);

    TTree tr1("tr1", "Training Sample");

    /*
     * Let's define training variables
     *
     * Units: Momentum in [GeV], angles in radians
     */

    double P_em1, th_em1, phi_em1;
    double P_em2, th_em2, phi_em2;
    double P_ep, th_ep, phi_ep;
    double Q2_em1, Q2_em2;

    tr1.Branch("P_em1", &P_em1, "P_em1/D");
    tr1.Branch("th_em1", &th_em1, "th_em1/D");
    tr1.Branch("phi_em1", &phi_em1, "phi_em1/D");
    tr1.Branch("P_em2", &P_em2, "P_em2/D");
    tr1.Branch("th_em2", &th_em2, "th_em2/D");
    tr1.Branch("phi_em2", &phi_em2, "phi_em2/D");
    tr1.Branch("P_ep", &P_ep, "P_ep/D");
    tr1.Branch("th_ep", &th_ep, "th_ep/D");
    tr1.Branch("phi_ep", &phi_ep, "phi_ep/D");
    tr1.Branch("Q2_em1", &Q2_em1, "Q2_em1/D");
    tr1.Branch("Q2_em2", &Q2_em2, "Q2_em2/D");

    constexpr int nMax_samePID = 20;
    constexpr int HTCC_TYPE = 15;
    constexpr double PCalEmin = 0.1;
    constexpr double PCalUmax = 405;
    constexpr double PCalUmin = 10;
    constexpr double PCalVmin = 10;
    constexpr double PCalWmin = 10;
    constexpr double vzMax = +2.;
    constexpr double vzMin = -8.;

    constexpr double vzMax_ep = 2;
    constexpr double vzMin_ep = -7;

    constexpr double vzMax_prot = 1.;
    constexpr double vzMin_prot = -6.;

    /*
     * -- The Beam electron index is 1
    * -- The Decay electron index is 3
    */

    constexpr int MC_ind_emBeam = 1;
    constexpr int MC_ind_emDecay = 3;
    constexpr int MC_ind_ep = 2;


    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;

    hipo::bank bRunConf(factory.getSchema("RUN::config"));
    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));
    hipo::bank bRecSC(factory.getSchema("REC::Scintillator"));
    hipo::bank bRecEV(factory.getSchema("REC::Event"));
    hipo::bank bMCPart(factory.getSchema("MC::Particle"));
    hipo::bank bMCRecMatch(factory.getSchema("MC::RecMatch"));

    int evCounter = 0;

    try {
        while(reader.next()) {

            evCounter = evCounter + 1;

            //if( evCounter > 155000 ){break;}
            if (evCounter % 10000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

            reader.read(event);
            event.getStructure(bRecPart);
            event.getStructure(bRecCalo);
            event.getStructure(bRecCC);
            event.getStructure(bRecSC);
            event.getStructure(bRecEV);
            event.getStructure(bRunConf);
            event.getStructure(bMCPart);

            /*
             * MC level dsitributions
             */
            double px_em_beam_MC = bMCPart.getFloat("px", MC_ind_emBeam);
            double py_em_beam_MC = bMCPart.getFloat("py", MC_ind_emBeam);
            double pz_em_beam_MC = bMCPart.getFloat("pz", MC_ind_emBeam);
            double p_em_beam_MC = sqrt(px_em_beam_MC*px_em_beam_MC + py_em_beam_MC*py_em_beam_MC + pz_em_beam_MC*pz_em_beam_MC);
            double th_em_beam_MC = acos(pz_em_beam_MC/p_em_beam_MC)*r2d;
            double phi_em_beam_MC = atan2(py_em_beam_MC, px_em_beam_MC)*r2d + 30.;
            phi_em_beam_MC = phi_em_beam_MC < 0 ? phi_em_beam_MC + 360: phi_em_beam_MC;

            double px_em_decay_MC = bMCPart.getFloat("px", MC_ind_emDecay);
            double py_em_decay_MC = bMCPart.getFloat("py", MC_ind_emDecay);
            double pz_em_decay_MC = bMCPart.getFloat("pz", MC_ind_emDecay);
            double p_em_decay_MC = sqrt(px_em_decay_MC*px_em_decay_MC + py_em_decay_MC*py_em_decay_MC + pz_em_decay_MC*pz_em_decay_MC);
            double th_em_decay_MC = acos(pz_em_decay_MC/p_em_decay_MC)*r2d;
            double phi_em_decay_MC = atan2(py_em_decay_MC, px_em_decay_MC)*r2d + 30.;
            phi_em_decay_MC = phi_em_decay_MC < 0 ? phi_em_decay_MC + 360: phi_em_decay_MC;

            double px_ep_MC = bMCPart.getFloat("px", MC_ind_ep);
            double py_ep_MC = bMCPart.getFloat("py", MC_ind_ep);
            double pz_ep_MC = bMCPart.getFloat("pz", MC_ind_ep);
            double p_ep_MC = sqrt(px_ep_MC*px_ep_MC + py_ep_MC*py_ep_MC + pz_ep_MC*pz_ep_MC);
            double th_ep_MC = acos(pz_ep_MC/p_ep_MC)*r2d;

            double Q2_beam_MC = 2*Eb*p_em_beam_MC*(1 - cos(th_em_beam_MC/r2d));
            double Q2_decay_MC = 2*Eb*p_em_decay_MC*(1 - cos(th_em_decay_MC/r2d));

            h_Q2_Decay_Beam1.Fill(Q2_beam_MC, Q2_decay_MC);

            h_th_p_MC_emBeam1.Fill(p_em_beam_MC, th_em_beam_MC);
            h_th_p_MC_emDecay1.Fill(p_em_decay_MC, th_em_decay_MC);
            h_th_p_MC_ep1.Fill(p_ep_MC, th_ep_MC);

            LVec L_em_decay_MC(px_em_decay_MC, py_em_decay_MC, pz_em_decay_MC, p_em_decay_MC);
            LVec L_em_beam_MC(px_em_beam_MC, py_em_beam_MC, pz_em_beam_MC, p_em_beam_MC);
            LVec L_ep_MC(px_ep_MC, py_ep_MC, pz_ep_MC, p_ep_MC);
            LVec L_emep_decay_MC = L_em_decay_MC + L_ep_MC;
            LVec L_emep_beam_MC  = L_em_beam_MC + L_ep_MC;

            double Minv_decay_MC = L_emep_decay_MC.M();
            double Minv_beam_MC = L_emep_beam_MC.M();

            h_Minv_Beam_Decay_MC1.Fill(Minv_decay_MC, Minv_beam_MC);
            h_P_em_Decay_Beam1.Fill(p_em_decay_MC, p_em_beam_MC);

            h_P_emep_Decay1.Fill(p_em_decay_MC, p_ep_MC);
            h_P_emep_Beam1.Fill(p_em_beam_MC, p_ep_MC);

            map<int, int> ind_HTCC;
            map<int, int> ind_PCal;
            map<int, int> ind_ECin;
            map<int, int> ind_ECout;

            int nPart = bRecPart.getRows();

            for ( int i_part = 0; i_part < nPart; i_part++ ) {
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

            int ind_em[nMax_samePID];
            int ind_ep[nMax_samePID];
            int ind_prot[nMax_samePID];
            int sec_em[nMax_samePID];
            int sec_ep[nMax_samePID];

            int n_em = 0;
            int n_ep = 0;
            int n_prot = 0;

            for ( int i_part = 0; i_part < nPart; i_part++ ) {
                RecParticle recp(bRecPart, bRecCalo, bRecCC, i_part, ind_PCal[i_part], ind_ECin[i_part], ind_ECout[i_part], ind_HTCC[i_part]);

                if (recp.pid() == 11 && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {

                    h_vz_em1.Fill(recp.vz());

                    h_PCal_lU_em1.Fill(recp.luPCal());
                    h_PCal_lV_em1.Fill(recp.lvPCal());
                    h_PCal_lW_em1.Fill(recp.lwPCal());

                    bool isPCalUMax = recp.luPCal() < PCalUmax;
                    bool isPCalVMin = recp.lvPCal() > PCalVmin;
                    bool isPCalWmin = recp.lwPCal() > PCalWmin;
                    bool isVZ = recp.vz() > vzMin && recp.vz() < vzMax;

                    if ( isPCalUMax && isPCalVMin && isPCalWmin && isVZ ) {
                        ind_em[n_em] = i_part;
                        n_em++;
                    }
                }else if ( recp.pid() == -11 && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000 ) {

                    h_vz_ep1.Fill(recp.vz());

                    bool isVZ = recp.vz() > vzMin_ep && recp.vz() < vzMax_ep;

                    h_PCal_lU_ep1.Fill(recp.luPCal());
                    h_PCal_lV_ep1.Fill(recp.lvPCal());
                    h_PCal_lW_ep1.Fill(recp.lwPCal());

                    bool isPCalUMin = recp.lvPCal() > PCalUmin;
                    bool isPCalVMin = recp.lvPCal() > PCalVmin;
                    bool isPCalWMin = recp.lwPCal() > PCalWmin;

                    if ( isVZ && isPCalUMin && isPCalVMin && isPCalWMin && ((recp.energyECin() < 0.001 && recp.SFPCal() > 0.11) || (recp.energyECin() >= 0.001 && (recp.SFPCal() + recp.SFECin()) > 0.17))  ) {
                        ind_ep[n_ep] = i_part;
                        n_ep++;
                    }
                } else if (recp.pid() == 2212) {

                    h_vz_prot.Fill(recp.vz());

                    bool isVZ = recp.vz() > vzMin_prot && recp.vz() < vzMax_prot;
                    if ( isVZ ) {
                        ind_prot[n_prot] = i_part;
                        n_prot++;
                    }
                }
            }

            if ( n_em == 2 && n_ep == 1 ) {

                /*
                 * em1(beam electron) and em2 will be assigned based on the MC::RecMatch bank
                 */

                event.getStructure(bMCRecMatch);
                int mc_ind_em1 = bMCRecMatch.getInt("mcindex", ind_em[0]);
                int mc_ind_em2 = bMCRecMatch.getInt("mcindex", ind_em[1]);
                double quality_em1 = bMCRecMatch.getFloat("quality", ind_em[0]);
                double quality_em2 = bMCRecMatch.getFloat("quality", ind_em[1]);

                h_quality_em2_em1_1.Fill(quality_em1, quality_em2);

                int ind_em_beam = -1;
                int ind_em_decay = -1;

                if ( mc_ind_em1 == MC_ind_emBeam ) {
                    ind_em_beam = ind_em[0];
                }else if ( mc_ind_em2 == MC_ind_emBeam) {
                    ind_em_beam = ind_em[1];
                }

                if ( mc_ind_em1 == MC_ind_emDecay ) {
                    ind_em_decay = ind_em[0];
                }else if ( mc_ind_em2 == MC_ind_emDecay) {
                    ind_em_decay = ind_em[1];
                }

                h_ind_embeam_decay1.Fill(ind_em_beam, ind_em_decay);

                if ( mc_ind_em1 == MC_ind_emBeam ) {
                    h_MatchQuality_em_beam1.Fill(quality_em1);
                    // if ( quality_em1 < 0.5 ) {
                    //     bMCRecMatch.show();
                    //     bRecPart.show();
                    //     bMCPart.show();
                    // }
                }else if ( mc_ind_em1 == MC_ind_emDecay ) {
                    h_MatchQuality_em_decay1.Fill(quality_em1);
                }


                RecParticle part_emBeam(bRecPart, bRecCalo, bRecCC, ind_em_beam, ind_PCal[ind_em_beam], ind_ECin[ind_em_beam], ind_ECout[ind_em_beam], ind_HTCC[ind_em_beam]);
                RecParticle part_emDecay(bRecPart, bRecCalo, bRecCC, ind_em_decay, ind_PCal[ind_em_decay], ind_ECin[ind_em_decay],ind_ECout[ind_em_decay],ind_HTCC[ind_em_decay]);
                RecParticle part_ep(bRecPart, bRecCalo, bRecCC, ind_ep[0], ind_PCal[ind_ep[0]], ind_ECin[ind_ep[0]],ind_ECout[ind_ep[0]],ind_HTCC[ind_ep[0]]);

                h_th_p_Rec_emDecay1.Fill( part_emDecay.p(), part_emDecay.th() );
                h_th_p_Rec_emBeam1.Fill( part_emBeam.p(), part_emBeam.th() );
                h_th_p_Rec_ep1.Fill( part_ep.p(), part_ep.th() );

                h_th_p_MC_emBeam2.Fill(p_em_beam_MC, th_em_beam_MC);
                h_th_p_MC_emDecay2.Fill(p_em_decay_MC, th_em_decay_MC);
                h_th_p_MC_ep2.Fill(p_ep_MC, th_ep_MC);

                h_Minv_Beam_Decay_MC2.Fill(Minv_decay_MC, Minv_beam_MC);
                h_Q2_Decay_Beam2.Fill(Q2_beam_MC, Q2_decay_MC);

                P_em1 = part_emBeam.p();
                th_em1 = part_emBeam.th()/r2d;
                phi_em1 = part_emBeam.phi()/r2d;
                P_em2 = part_emDecay.p();
                th_em2 = part_emDecay.th()/r2d;
                phi_em2 = part_emDecay.phi()/r2d;
                P_ep = part_ep.p();
                th_ep = part_ep.th()/r2d;
                phi_ep = part_ep.phi()/r2d;

                Q2_em1 = 2*Eb*P_em1*(1 - cos(th_em1));
                Q2_em2 = 2*Eb*P_em2*(1 - cos(th_em2));

                if ( ind_em_beam>= 0 && ind_em_decay >= 0 ) {

                    LVec L_em1(part_emBeam.px(), part_emBeam.py(), part_emBeam.pz(), part_emBeam.p());
                    LVec L_em2(part_emDecay.px(), part_emDecay.py(), part_emDecay.pz(), part_emDecay.p());
                    LVec L_ep(part_ep.px(), part_ep.py(), part_ep.pz(), part_ep.p());
                    ddvcs_kine_Noprot.SetKineem1em2ep(L_em1, L_em2, L_ep);

                    h_P_beam_em_Rec_MC1.Fill(p_em_beam_MC, part_emBeam.p());
                    h_th_beam_em_Rec_MC1.Fill(th_em_beam_MC, part_emBeam.th());
                    h_phi_beam_em_Rec_MC1.Fill(phi_em_beam_MC, part_emBeam.phi());
                    h_P_decay_em_Rec_MC1.Fill(p_em_decay_MC, part_emDecay.p());
                    h_th_decay_em_Rec_MC1.Fill(th_em_decay_MC, part_emDecay.th());
                    h_phi_decay_em_Rec_MC1.Fill(phi_em_decay_MC, part_emDecay.phi());

                    double Mx2 = ddvcs_kine_Noprot.GetMxRecoil();
                    h_Mmis1.Fill(Mx2);

                    if ( Mx2 > MxRecoil_Min && Mx2 < MxRecoil_Max ) {
                        tr1.Fill();
                    }
                }
            }

        }
    }catch(exception& e) {
        cerr << "Exception thrown" <<e.what()<< endl;
    }

    gDirectory->Write();
    outFile.Close();
    return EXIT_SUCCESS;
}