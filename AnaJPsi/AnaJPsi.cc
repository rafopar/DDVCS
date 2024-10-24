/* 
 * File:   AnaJPsi.cc
 * Author: rafopar
 *
 * Created on October 21, 2024, 5:25 PM
 */


#include <ctime>
#include <chrono>
#include <cstdlib>

#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include <RecParticle.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>

using namespace std;

void InitFunction(TF1*, std::string kewWord);

/*
 * Rotates the LorentzVector L by an angle "angle", n a way
 * that the azimuthal angle is unchanged.
 * angle is in degrees
 */
void RotateTheta(TLorentzVector& L, double angle);

struct MCPart {
    int pid;
    double px;
    double py;
    double pz;
    double p;
    double vx;
    double vy;
    double vz;
};

/*
 * In this analysis we will check how wide will be the J/psi mass after all corrections 
 * we apply with the uCLAS12 setup.
 */
int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Please provide the run number." << endl;
        cout << "Exiting" << endl;
    }

    const int run = atoi(argv[1]);
    const int HTCC_TYPE = 15;
    const double Mmu = 0.10566;
    const double Mp = 0.9383;

    const double th_MuMax = 40; // We ignore angles above this value
    
    std::map<int, double> m_Eb;

    m_Eb[1] = 10.6;
    m_Eb[2] = 22;

    double Eb = m_Eb[run];


    /*
     * Th is functions looks like works well for Eloss: [0] + [1]*x + [2]/x + [3]/(x*x) + [4]/(x*x*x)
     */
    TF1 *f_Eloss = new TF1("f_Eloss", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 80.);
    f_Eloss->SetNpx(4500);

    TF1 *f_ThCorr = new TF1("f_ThCorr", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0])) + [5]/((x-[0])*(x-[0])*(x-[0]))", 0., 80.);
    f_ThCorr->SetNpx(4500);

    TF1 *f_PhiCorr = new TF1("f_PhiCorr", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 80.);
    f_PhiCorr->SetNpx(4500);

    TF1 *f_Elos_mup = (TF1*) f_Eloss->Clone("f_Elos_mup");
    TF1 *f_Elos_mum = (TF1*) f_Eloss->Clone("f_Elos_mum");
    TF1 *f_ThCorr_mup = (TF1*) f_ThCorr->Clone("f_ThCorr_mup");
    TF1 *f_ThCorr_mum = (TF1*) f_ThCorr->Clone("f_ThCorr_mum");
    TF1 *f_PhiCorr_mup = (TF1*) f_PhiCorr->Clone("f_PhiCorr_mup");
    TF1 *f_PhiCorr_mum = (TF1*) f_PhiCorr->Clone("f_PhiCorr_mum");

    /*
     * Now let's initialize these functions with paramaters
     */

    InitFunction(f_Elos_mup, Form("mup_ElossFunc_Pars_Run_%d", run));
    InitFunction(f_Elos_mum, Form("mum_ElossFunc_Pars_Run_%d", run));
    InitFunction(f_ThCorr_mup, Form("mup_DeltaTheta_Corr_Run_%d", run));
    InitFunction(f_ThCorr_mum, Form("mum_DeltaTheta_Corr_Run_%d", run));
    InitFunction(f_PhiCorr_mup, Form("mup_DeltaPhi_Corr_Run_%d", run));
    InitFunction(f_PhiCorr_mum, Form("mum_DeltaPhi_Corr_Run_%d", run));

    TLorentzVector L_em_mc, L_mum_mc, L_mup_mc, L_p_mc;
    TLorentzVector L_em_rec, L_mum_rec, L_mup_rec; // Recon LorentzVectors
    TLorentzVector L_em_cor, L_mum_cor, L_mup_cor; // Corrected LorentzVectors (Elos, delta_theta and delta_phi)
    TLorentzVector L_beam, L_targ;
    L_beam.SetPxPyPzE(0., 0., Eb, Eb);
    L_targ.SetPxPyPzE(0., 0., 0., Mp);

    TFile file_out(Form("AnaJPsi_Run_%d.root", run), "Recreate");

    TH2D h_th_P_mup1("h_th_P_mup1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mup2("h_th_P_mup2", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mum1("h_th_P_mum1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mum2("h_th_P_mum2", "", 200, 0., 1.05 * Eb, 200, 0., 45);

    TH2D h_dPP_P_mup1("h_dPP_P_mup1", "", 200, 0., 1.05 * Eb, 200, -0.9, 0.05);
    TH2D h_dPP_P_mum1("h_dPP_P_mum1", "", 200, 0., 1.05 * Eb, 200, -0.9, 0.05);

    TH2D h_DeltaP_P_mup1("h_DeltaP_P_mup1", "", 200, 0., 1.05 * Eb, 200, -2, 0.05);
    TH2D h_DeltaP_P_mum1("h_DeltaP_P_mum1", "", 200, 0., 1.05 * Eb, 200, -2, 0.05);

    TH2D h_DeltaTheta_P_mup1("h_DeltaTheta_P_mup1", "", 200, 0., 1.05 * Eb, 200, -20., 20);
    TH2D h_DeltaTheta_P_mum1("h_DeltaTheta_P_mum1", "", 200, 0., 1.05 * Eb, 200, -20., 20);

    TH2D h_DeltaTheta_Theta_mup1("h_DeltaTheta_Theta_mup1", "", 200, 0., 65, 200, -20., 40);
    TH2D h_DeltaTheta_Theta_mum1("h_DeltaTheta_Theta_mum1", "", 200, 0., 65, 200, -20., 40);

    TH2D h_DeltaPhi_Pt_mup1("h_DeltaPhi_Pt_mup1", "", 200., 0., 0.18*Eb, 200, -180., 180.);
    TH2D h_DeltaPhi_Pt_mum1("h_DeltaPhi_Pt_mum1", "", 200., 0., 0.18*Eb, 200, -180., 180.);

    TH1D h_Mmumu_Rec1("h_Mmumu_Rec1", "", 200, 1., 4.);
    TH1D h_Mmumu_Rec_Corr1("h_Mmumu_Rec_Corr1", "", 200, 1., 4.);
    TH1D h_Mmumu_Rec_Corr2("h_Mmumu_Rec_Corr2", "", 200, 1., 4.);
    TH1D h_Mmumu_Rec_Corr3("h_Mmumu_Rec_Corr3", "", 200, 1., 4.);
    TH1D h_Mmumu_Rec_AngleFixCorr1("h_Mmumu_Rec_AngleFixCorr1", "", 200, 1., 4.);

    TH2D h_th_P_mup_MC1("h_th_P_mup_MC1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mup_MC2("h_th_P_mup_MC2", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mum_MC1("h_th_P_mum_MC1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mum_MC2("h_th_P_mum_MC2", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    
    
    hipo::reader reader;
    //reader.open("Data/Skim_ZeroSuppr_2247_All.hipo");
    reader.open(Form("Data/JPsi_Run%d_All.hipo", run));

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

    // Start time
    auto start = std::chrono::high_resolution_clock::now();
    auto startTime = std::chrono::system_clock::to_time_t(start);
    cout << "The loop starts at " << std::ctime(&startTime) << endl;

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


            int n_mum = 0;
            int n_mup = 0;
            int n_prot = 0;
            int n_charged = 0;

            vector<RecParticle> v_recp_mum;
            vector<RecParticle> v_recp_mup;
            vector<RecParticle> v_recp_p;


            for (int ipart = 0; ipart < nPart; ipart++) {
                RecParticle recp(bRecPart, bRecCalo, bRecCC, ipart, ind_PCal[ipart], ind_ECin[ipart], ind_ECout[ipart], ind_HTCC[ipart]);

                if ( (recp.pid() == 211 || recp.pid() == -13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000 && recp.th() < th_MuMax ) {
                    v_recp_mup.push_back(recp);
                    h_th_P_mup1.Fill(recp.p(), recp.th());
                    n_mup = n_mup + 1;
                } else if ( (recp.pid() == -211 || recp.pid() == 13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000 && recp.th() < th_MuMax ) {
                    v_recp_mum.push_back(recp);
                    h_th_P_mum1.Fill(recp.p(), recp.th());
                    n_mum = n_mum + 1;
                }
            }

            MCPart mc_mum, mc_mup, mc_p;
            for (int imc = 0; imc < nMCPart; imc++) {
                MCPart curPart;

                /*
                 *  For some reason the MC::Particle bank has several more entries w/ pid == 11 and p = 11 GeV
                 *   It has px and py = 0, so using that, we can skip those particles.
                 */
                if (bMCPart.getFloat("px", imc) == 0 && bMCPart.getFloat("py", imc) == 0) {
                    continue;
                }

                curPart.pid = bMCPart.getInt("pid", imc);
                curPart.px = double(bMCPart.getFloat("px", imc));
                curPart.py = double(bMCPart.getFloat("py", imc));
                curPart.pz = double(bMCPart.getFloat("pz", imc));
                curPart.p = sqrt(curPart.px * curPart.px + curPart.py * curPart.py + curPart.pz * curPart.pz);
                curPart.vx = double(bMCPart.getFloat("vx", imc));
                curPart.vy = double(bMCPart.getFloat("vy", imc));
                curPart.vz = double(bMCPart.getFloat("vz", imc));

                if (curPart.pid == 13) {
                    mc_mum = curPart;
                } else if (curPart.pid == -13) {
                    mc_mup = curPart;
                } else if (curPart.pid == 2212) {
                    mc_p = curPart;
                }
            }

            L_mum_mc.SetPxPyPzE(mc_mum.px, mc_mum.py, mc_mum.pz, sqrt(mc_mum.p * mc_mum.p + Mmu * Mmu));
            
            L_mup_mc.SetPxPyPzE(mc_mup.px, mc_mup.py, mc_mup.pz, sqrt(mc_mup.p * mc_mup.p + Mmu * Mmu));

            h_th_P_mup_MC1.Fill( L_mup_mc.P(), L_mup_mc.Theta()*TMath::RadToDeg() );
            h_th_P_mum_MC1.Fill( L_mum_mc.P(), L_mum_mc.Theta()*TMath::RadToDeg() );

            if (n_mum == 1 && n_mup == 1) {

                RecParticle rec_mup = v_recp_mup.at(0);
                RecParticle rec_mum = v_recp_mum.at(0);

                double th_mup_Rec = rec_mup.th();
                double th_mum_Rec = rec_mum.th();

                double phi_mup_Rec = rec_mup.phi() * TMath::RadToDeg();
                double phi_mum_Rec = rec_mum.phi() * TMath::RadToDeg();

                double th_mup_MC = L_mup_mc.Theta() * TMath::RadToDeg();
                double th_mum_MC = L_mum_mc.Theta() * TMath::RadToDeg();

                double pt_mup_Rec = sqrt(rec_mup.px() * rec_mup.px() + rec_mup.py() * rec_mup.py());
                double pt_mum_Rec = sqrt(rec_mum.px() * rec_mum.px() + rec_mum.py() * rec_mum.py());

                double phi_mup_rec = rec_mup.phi() - 30;
                phi_mup_rec = phi_mup_rec < 0 ? phi_mup_rec + 360 : phi_mup_rec;
                double phi_mum_rec = rec_mum.phi() - 30;
                phi_mum_rec = phi_mum_rec < 0 ? phi_mum_rec + 360 : phi_mum_rec;

                double deltaPhi_mup = phi_mup_rec - L_mup_mc.Phi() * TMath::RadToDeg();
                deltaPhi_mup = TMath::Abs(deltaPhi_mup) < 180 ? deltaPhi_mup : deltaPhi_mup - 360 * (TMath::Abs(deltaPhi_mup) / deltaPhi_mup);

                double deltaPhi_mum = phi_mum_rec - L_mum_mc.Phi() * TMath::RadToDeg();
                deltaPhi_mum = TMath::Abs(deltaPhi_mum) < 180 ? deltaPhi_mum : deltaPhi_mum - 360 * (TMath::Abs(deltaPhi_mum) / deltaPhi_mum);

                h_th_P_mup2.Fill(rec_mup.p(), th_mup_Rec);
                h_th_P_mum2.Fill(rec_mum.p(), th_mum_Rec);

                h_dPP_P_mup1.Fill(rec_mup.p(), (rec_mup.p() - L_mup_mc.P()) / L_mup_mc.P());
                h_dPP_P_mum1.Fill(rec_mum.p(), (rec_mum.p() - L_mum_mc.P()) / L_mum_mc.P());

                h_DeltaP_P_mup1.Fill(rec_mup.p(), rec_mup.p() - L_mup_mc.P());
                h_DeltaP_P_mum1.Fill(rec_mum.p(), rec_mum.p() - L_mum_mc.P());

                h_DeltaTheta_P_mup1.Fill(rec_mup.p(), th_mup_Rec - th_mup_MC);
                h_DeltaTheta_P_mum1.Fill(rec_mum.p(), th_mum_Rec - th_mum_MC);

                h_DeltaTheta_Theta_mup1.Fill(th_mup_Rec, th_mup_Rec - th_mup_MC);
                h_DeltaTheta_Theta_mum1.Fill(th_mum_Rec, th_mum_Rec - th_mum_MC);

                h_DeltaPhi_Pt_mup1.Fill(pt_mup_Rec, deltaPhi_mup);
                h_DeltaPhi_Pt_mum1.Fill(pt_mum_Rec, deltaPhi_mum);

                L_mup_rec.SetPxPyPzE(rec_mup.px(), rec_mup.py(), rec_mup.pz(), sqrt(rec_mup.p() * rec_mup.p() + Mmu * Mmu));
                L_mum_rec.SetPxPyPzE(rec_mum.px(), rec_mum.py(), rec_mum.pz(), sqrt(rec_mum.p() * rec_mum.p() + Mmu * Mmu));

                h_th_P_mup_MC2.Fill(L_mup_mc.P(), L_mup_mc.Theta() * TMath::RadToDeg());
                h_th_P_mum_MC2.Fill(L_mum_mc.P(), L_mum_mc.Theta() * TMath::RadToDeg());

                TLorentzVector L_mumu = L_mum_rec + L_mup_rec;
                
                double m_mumu = L_mumu.M();
                h_Mmumu_Rec1.Fill(m_mumu);
                
                /*
                 * Let's make kinematic corrections below.
                 */

                double p_mup_corr = L_mup_rec.P() / (1 + f_Elos_mup->Eval(L_mup_rec.P()));
                double momScale = p_mup_corr / L_mup_rec.P();
                L_mup_cor.SetPxPyPzE(L_mup_rec.Px() * momScale, L_mup_rec.Py() * momScale, L_mup_rec.Pz() * momScale, sqrt(p_mup_corr * p_mup_corr + Mmu * Mmu));

                // LorentzVector that has only Momentum corrected, but angles are fixed to generated values
                TLorentzVector L_mup_AngleFixCor = L_mup_cor;
                L_mup_AngleFixCor.SetTheta( L_mup_mc.Theta() );
                L_mup_AngleFixCor.SetPhi( L_mup_mc.Phi() );
                
                double p_mum_corr = L_mum_rec.P() / (1 + f_Elos_mum->Eval(L_mum_rec.P()));
                momScale = p_mum_corr / L_mum_rec.P();
                L_mum_cor.SetPxPyPzE(L_mum_rec.Px() * momScale, L_mum_rec.Py() * momScale, L_mum_rec.Pz() * momScale, sqrt(p_mum_corr * p_mum_corr + Mmu * Mmu));

                // LorentzVector that has only Momentum corrected, but angles are fixed to generated values
                TLorentzVector L_mum_AngleFixCor = L_mum_cor;
                L_mum_AngleFixCor.SetTheta( L_mum_mc.Theta() );
                L_mum_AngleFixCor.SetPhi( L_mum_mc.Phi() );
                
                TLorentzVector L_mumu_corr1 = L_mum_cor + L_mup_cor;
                double m_mumu_corr1 = L_mumu_corr1.M();
                h_Mmumu_Rec_Corr1.Fill(m_mumu_corr1);
                
                
                double mup_delta_th_corr = f_ThCorr_mup->Eval(th_mup_Rec);
                RotateTheta(L_mup_cor, mup_delta_th_corr);

                double mum_delta_th_corr = f_ThCorr_mum->Eval(th_mum_Rec);
                RotateTheta(L_mum_cor, mum_delta_th_corr);
                
                TLorentzVector L_mumu_corr2 = L_mum_cor + L_mup_cor;
                double m_mumu_corr2 = L_mumu_corr2.M();
                h_Mmumu_Rec_Corr2.Fill(m_mumu_corr2);
                
                double mup_delta_Phi = -f_PhiCorr_mup->Eval(pt_mup_Rec);               
                L_mup_cor.RotateZ(mup_delta_Phi * TMath::DegToRad());

                double mum_delta_Phi = -f_PhiCorr_mum->Eval(pt_mum_Rec);
                L_mum_cor.RotateZ(mum_delta_Phi * TMath::DegToRad());
                
                TLorentzVector L_mumu_corr3 = L_mum_cor + L_mup_cor;
                double m_mumu_corr3 = L_mumu_corr3.M();
                h_Mmumu_Rec_Corr3.Fill(m_mumu_corr3);
                
                TLorentzVector L_mumu_AngleFixCorr = L_mup_AngleFixCor + L_mum_AngleFixCor;
                double m_mumu_AngleFixXorr = L_mumu_AngleFixCorr.M();
                
                h_Mmumu_Rec_AngleFixCorr1.Fill(m_mumu_AngleFixXorr);
            }

        }
    } catch (const char msg) {
        cerr << msg << endl;
    }

    gDirectory->Write();

    return 0;
}

void InitFunction(TF1* f, std::string keyword) {

    std::string fname = Form("Pars/%s.dat", keyword.c_str());

    ifstream inp_dat(fname);

    if (!inp_dat.is_open()) {
        cout << "Can not open the file " << fname << ". Exiting" << endl;
        exit(1);
    }

    while (!inp_dat.eof()) {
        int ind;
        double par;

        inp_dat>>ind;
        inp_dat>>par;
        f->SetParameter(ind, par);
    }
}

void RotateTheta(TLorentzVector& L, double theta) {

    TVector3 momentum_proj_xy(L.Px(), L.Py(), 0.0);

    TVector3 rotation_axis = momentum_proj_xy.Cross(L.Vect());

    L.Rotate(theta * TMath::DegToRad(), rotation_axis);
}