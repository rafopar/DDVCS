/* 
 * File:   AnaGEMC.cc
 * Author: rafopar
 *
 * Created on August 30, 2024, 1:08 PM
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

/*
 *  In this analysis, we should reconstruct muon pairs from the REC::Particle bank,
 *  while for electrons we will still apply geometric acceptance cuts.
 * 
 * ************** What needs to be studied *******************
 *   * Acceptance of both muons
 *   * How much energy muons lose passing through the Calorimeter and tungsten shielding
 *   * How good (or better to say bad) is vertex reconstruction
 *   * How close are reconstructed angles to generated angles
 */

bool emAcc(TLorentzVector&);

void InitFunction(TF1*, std::string kewWord);

/*
 * Rotates the LorentzVector L by an angle "angle", n a way
 * that the azimuthal angle is unchanged.
 * angle is in degrees
 */
void RotateTheta(TLorentzVector& L, double angle);

int findX_Xi_Bin(double x, double xi);

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
 * 
 */
int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Please provide the run number." << endl;
        cout << "Exiting" << endl;
    }

    const int HTCC_TYPE = 15;
    const double Mmu = 0.10566;
    const double Mp = 0.9383;
    const double r2d = TMath::RadToDeg();

    const double InstLumin = 1.e37;
    const double secondPerDay = 3600 * 24;
    const double nDays = 200;
    const double totLumi = InstLumin * secondPerDay*nDays;
    const double pbn = 1.e-36;

    const double MM2_max = 1.5;
    const double MM2_min = 0.4;

    const double bin0_x_max = -0.07;
    const double bin0_x_min = -0.12;
    const double bin0_xi_max = 0.24;
    const double bin0_xi_min = 0.19;

    const double bin1_x_max = 0.08;
    const double bin1_x_min = 0.02;
    const double bin1_xi_max = 0.24;
    const double bin1_xi_min = 0.18;

    const double th_MuMin_MC = 7.1; // We ignore angles below this value
    const double th_MuMax = 40; // We ignore angles above this value

    const int nMaxPCalStripsPos = 5;
    const int nMaxECinStripsPos = 5;
    const int nMaxECoutStripsPos = 5;
    const int nMaxPCalStripsNeg = 7;
    const int nMaxECinStripsNeg = 5;
    const int nMaxECoutStripsNeg = 5;


    const int n_xi_x_bins = 2;
    const int n_xi_x_Q2bins = 2;

    int run = atoi(argv[1]);

    std::map<int, double> m_Eb;
    m_Eb[17] = 10.6;
    m_Eb[1701] = 10.6;
    m_Eb[18] = 22.;
    m_Eb[19] = 10.6;
    m_Eb[20] = 22.;
    m_Eb[9113] = 10.6;
    m_Eb[9114] = 10.6;
    m_Eb[9115] = 22.;
    m_Eb[9116] = 22.;
    m_Eb[1022] = 10.6; // This is a merged run of Run1 Inclusive AND Run 23 from GRAPE
    m_Eb[1023] = 10.6; // This is a merged run of Run1 Inclusive AND Run 23 from GRAPE
    m_Eb[2022] = 10.6; // This is a merged run of Run1 Inclusive AND Run 23 from GRAPE
    m_Eb[2023] = 10.6; // This is a merged run of Run1 Inclusive AND Run 23 from GRAPE

    std::map<int, int> m_MC_mupPID;
    std::map<int, int> m_MC_mumPID;
    m_MC_mupPID[9114] = 211;
    m_MC_mumPID[9114] = -211;
    m_MC_mupPID[9116] = 211;
    m_MC_mumPID[9116] = -211;
    m_MC_mupPID[9113] = -13;
    m_MC_mumPID[9113] = 13;
    m_MC_mupPID[9115] = -13;
    m_MC_mumPID[9115] = 13;
    m_MC_mupPID[17] = -13;
    m_MC_mumPID[17] = 13;
    m_MC_mupPID[1701] = 211;
    m_MC_mumPID[1701] = -211;
    m_MC_mupPID[18] = -13;
    m_MC_mumPID[18] = 13;
    m_MC_mupPID[19] = -13;
    m_MC_mumPID[19] = 13;
    m_MC_mupPID[20] = -13;
    m_MC_mumPID[20] = 13;
    m_MC_mupPID[1022] = -13;
    m_MC_mumPID[1022] = 13;
    m_MC_mupPID[1023] = -13;
    m_MC_mumPID[1023] = 13;
    m_MC_mupPID[2022] = -13;
    m_MC_mumPID[2022] = 13;
    m_MC_mupPID[2023] = -13;
    m_MC_mumPID[2023] = 13;

    std::map<int, double> m_tot_xSec;
    m_tot_xSec[17] = 0.0642577; // in pB
    m_tot_xSec[1701] = 0.0642577; // in pB
    m_tot_xSec[18] = 0.0816628; // in pB
    m_tot_xSec[19] = 0.0104044; // in pB
    m_tot_xSec[20] = 0.0263213; // in pB
    m_tot_xSec[9113] = 1; // in pB, for Rho, the number is from Harut
    m_tot_xSec[9114] = 10000.; // in pB, This is Rho with pi-pi+
    m_tot_xSec[9115] = 0.9; // in pB, This is Rho with mu-mu+
    m_tot_xSec[9116] = 9000.; // in pB, This is Rho with mu-mu+
    m_tot_xSec[1022] = 0.0026611735; // in pB, Obtained as sigma_1*sigma_2*4ns*10^{37}
    m_tot_xSec[1023] = 0.016686376; // in pB, Obtained as sigma_1*sigma_2*4ns*10^{37}
    m_tot_xSec[2022] = 0.0042309827; // in pB, Obtained as sigma_1*sigma_2*4ns*10^{37}
    m_tot_xSec[2023] = 0.0265295262; // in pB, Obtained as sigma_1*sigma_2*4ns*10^{37}

    std::map<int, double> m_Qp2_Q2Cut_Offset_bin0;
    std::map<int, double> m_Qp2_Q2Cut_Slope_bin0;
    std::map<int, double> m_Qp2_Q2Cut_Offset_bin1;
    std::map<int, double> m_Qp2_Q2Cut_Slope_bin1;
    // ****** 10.6 GeV *******
    m_Qp2_Q2Cut_Offset_bin0[17] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[1701] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[19] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[9113] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[9114] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[1023] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[1022] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[2023] = 4.6;
    m_Qp2_Q2Cut_Offset_bin0[2022] = 4.6;
    m_Qp2_Q2Cut_Slope_bin0[17] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[19] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[1022] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[1023] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[2022] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[2023] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[9113] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[9114] = -1. / 1.65;


    m_Qp2_Q2Cut_Offset_bin1[17] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[1701] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[19] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[9113] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[9114] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[1022] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[1023] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[2022] = 6.9;
    m_Qp2_Q2Cut_Offset_bin1[2023] = 6.9;
    m_Qp2_Q2Cut_Slope_bin1[17] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[19] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[9113] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[9114] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[1022] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[1023] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[2022] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[2023] = -1. / 0.55;
    // ****** 22 GeV *******
    m_Qp2_Q2Cut_Offset_bin0[18] = 11.;
    m_Qp2_Q2Cut_Offset_bin0[20] = 11.;
    m_Qp2_Q2Cut_Offset_bin0[9115] = 11.;
    m_Qp2_Q2Cut_Offset_bin0[9116] = 11.;
    m_Qp2_Q2Cut_Slope_bin0[18] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[20] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[9115] = -1. / 1.65;
    m_Qp2_Q2Cut_Slope_bin0[9116] = -1. / 1.65;


    m_Qp2_Q2Cut_Offset_bin1[18] = 13.5;
    m_Qp2_Q2Cut_Offset_bin1[20] = 13.5;
    m_Qp2_Q2Cut_Offset_bin1[9115] = 13.5;
    m_Qp2_Q2Cut_Offset_bin1[9116] = 13.5;
    m_Qp2_Q2Cut_Slope_bin1[18] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[20] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[9115] = -1. / 0.55;
    m_Qp2_Q2Cut_Slope_bin1[9116] = -1. / 0.55;

    double Qp2Q2_bin0_Offset = m_Qp2_Q2Cut_Offset_bin0[run];
    double Qp2Q2_bin0_Slope = m_Qp2_Q2Cut_Slope_bin0[run];
    double Qp2Q2_bin1_Offset = m_Qp2_Q2Cut_Offset_bin1[run];
    double Qp2Q2_bin1_Slope = m_Qp2_Q2Cut_Slope_bin1[run];


    /*
     * For the t Cut we use 0.5 GeV2 for 12 GeV, and 1 GeV for 22 GeV
     */
    std::map<int, double> m_tMCut;
    m_tMCut[17] = 0.5;
    m_tMCut[19] = 0.5;
    m_tMCut[9113] = 0.5;
    m_tMCut[9114] = 0.5;
    m_tMCut[18] = 1.;
    m_tMCut[20] = 1.;
    m_tMCut[9115] = 1.;
    m_tMCut[9116] = 1.;

    const double t_Max_xi_xiStudy = m_tMCut[run];


    const double MinvMin1 = 1.2;
    const double MinvMin2 = 1.8;

    if (m_tot_xSec.find(run) == m_tot_xSec.end() || m_Eb.find(run) == m_Eb.end()) {
        cout << "The xsec or the beam energy is not set for the run " << run << "." << endl;
        cout << "Exiting." << endl;
        exit(1);
    }

    /*
     * Most of the time muon id should be 13/-13, but for some studies
     * we need pions. So by default we will assume 13/-13
     */
    int MC_mupPID = m_MC_mupPID[run];
    int MC_mumPID = m_MC_mumPID[run];

    double tot_xSec = m_tot_xSec[run];

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

    /*
     * Function(s) below are for checking quick things, just test functions
     */

    TF1 *f_Pol4Test1 = new TF1("f_Pol4Test1", "[0] + x*([1] + x*( [2] + x*([3] + x*[4]) ) )", 0., 22.);
    f_Pol4Test1->SetParameters(-0.781752, 0.118757, -0.00456906, -0.000312771, 1.88916e-05);
    TF1 *f_Pol3Test2 = new TF1("f_Pol3Test2", "pol3", 0., 22.);
    f_Pol3Test2->SetParameters(-0.384326, 0.059916, -0.00407646, 0.000100553);

    TLorentzVector L_em_mc, L_mum_mc, L_mup_mc, L_p_mc;
    TLorentzVector L_em_rec, L_mum_rec, L_mup_rec; // Recon LorentzVectors
    TLorentzVector L_em_cor, L_mum_cor, L_mup_cor; // Corrected LorentzVectors (Elos, delta_theta and delta_phi)
    TLorentzVector L_beam, L_targ;
    L_beam.SetPxPyPzE(0., 0., Eb, Eb);
    L_targ.SetPxPyPzE(0., 0., 0., Mp);

    TFile file_out(Form("AnaDDVCS_Run_%d.root", run), "Recreate");

    TH2D h_th_P_mup1("h_th_P_mup1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mup2("h_th_P_mup2", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mum1("h_th_P_mum1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mum2("h_th_P_mum2", "", 200, 0., 1.05 * Eb, 200, 0., 45);


    TH2D h_MC_th_phi_mup_miss1("h_MC_th_phi_mup_miss1", "", 200, -180., 180., 200, 0., 60.);
    TH2D h_MC_th_phi_mup_det1("h_MC_th_phi_mup_det1", "", 200, -180., 180., 200, 0., 60.);
    TH2D h_MC_th_phi_mum_miss1("h_MC_th_phi_mum_miss1", "", 200, -180., 180., 200, 0., 60.);
    TH2D h_MC_th_phi_mum_det1("h_MC_th_phi_mum_det1", "", 200, -180., 180., 200, 0., 60.);

    TH2D h_MC_th_P_mup_miss1("h_MC_th_P_mup_miss1", "", 200, 0., 1.05 * Eb, 200, 0., 45.);
    TH2D h_MC_th_P_mup_det1("h_MC_th_P_mup_det1", "", 200, 0., 1.05 * Eb, 200, 0., 45.);
    TH2D h_MC_th_P_mum_miss1("h_MC_th_P_mum_miss1", "", 200, 0., 1.05 * Eb, 200, 0., 45.);
    TH2D h_MC_th_P_mum_det1("h_MC_th_P_mum_det1", "", 200, 0., 1.05 * Eb, 200, 0., 45.);
    TH2D h_MC_th_P_em_miss1("h_MC_th_P_em_miss1", "", 200, 0., 1.05 * Eb, 200, 0., 45.);
    TH2D h_MC_th_P_em_det1("h_MC_th_P_em_det1", "", 200, 0., 1.05 * Eb, 200, 0., 45.);
    
    TH2D h_MC_th_P_em_FS1("h_MC_th_P_em_FS1", "", 200, 0., 1.05 * Eb, 200, 0., 45.); // In the e-mu-mu+ Final state
    TH2D h_MC_th_P_em_FS2("h_MC_th_P_em_FS2", "", 200, 0., 1.05 * Eb, 200, 0., 45.); // In the e-mu-mu+ Final state, after Mx2 cut
    TH2D h_MC_th_P_em_FS3("h_MC_th_P_em_FS3", "", 200, 0., 1.05 * Eb, 200, 0., 45.); // In the e-mu-mu+ Final state, after the Mx2 and tM cut

    TH2D h_MC_th_P_em1("h_MC_th_P_em1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_MC_th_P_em2("h_MC_th_P_em2", "", 200, 0., 1.05 * Eb, 200, 0., 45);

    TH2D h_dPP_P_mup1("h_dPP_P_mup1", "", 200, 0., 1.05 * Eb, 200, -0.9, 0.05);
    TH2D h_dPP_P_mum1("h_dPP_P_mum1", "", 200, 0., 1.05 * Eb, 200, -0.9, 0.05);

    /*
     * Those with a suffix Test indicate some diagnostic histograms which will be used
     * to check some staff in a short term
     */
    TH2D h_th_P_mupTEST1("h_th_P_mupTEST1", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mupTEST2("h_th_P_mupTEST2", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mupTEST3("h_th_P_mupTEST3", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_th_P_mupTEST4("h_th_P_mupTEST4", "", 200, 0., 1.05 * Eb, 200, 0., 45);
    TH2D h_dPP_P_mup_Test1("h_dPP_P_mup_Test1", "", 200, 0., 1.05 * Eb, 200, -0.9, 0.05);
    TH2D h_dPP_P_mup_Test2("h_dPP_P_mup_Test2", "", 200, 0., 1.05 * Eb, 200, -0.9, 0.05);


    TH2D h_DeltaP_P_mup1("h_DeltaP_P_mup1", "", 200, 0., 1.05 * Eb, 200, -2, 0.05);
    TH2D h_DeltaP_P_mum1("h_DeltaP_P_mum1", "", 200, 0., 1.05 * Eb, 200, -2, 0.05);

    TH2D h_DeltaTheta_P_mup1("h_DeltaTheta_P_mup1", "", 200, 0., 1.05 * Eb, 200, -20., 20);
    TH2D h_DeltaTheta_P_mum1("h_DeltaTheta_P_mum1", "", 200, 0., 1.05 * Eb, 200, -20., 20);

    TH2D h_DeltaTheta_Theta_mup1("h_DeltaTheta_Theta_mup1", "", 200, 0., 65, 200, -20., 40);
    TH2D h_DeltaTheta_Theta_mum1("h_DeltaTheta_Theta_mum1", "", 200, 0., 65, 200, -20., 40);

    TH2D h_DeltaPhi_Pt_mup1("h_DeltaPhi_Pt_mup1", "", 200., 0., 1.8, 200, -180., 180.);
    TH2D h_DeltaPhi_Pt_mum1("h_DeltaPhi_Pt_mum1", "", 200., 0., 1.8, 200, -180., 180.);

    TH1D h_Mmis1("h_Mmis1", "", 200, -1., 6.);

    TH1D h_Mmis_corr1("h_Mmis_corr1", "", 200, -1., 6.);
    TH1D h_Mmis_corr2("h_Mmis_corr2", "", 200, -1., 6.);
    TH1D h_Mmis_corr3("h_Mmis_corr3", "", 200, -1., 6.);
    TH1D h_Mmis_AngleFixcorr_WideRange1("h_Mmis_AngleFixcorr_WideRange1", "", 200, -14., 6.);
    TH1D h_Mmis_AngleFixcorr1("h_Mmis_AngleFixcorr1", "", 200, -1., 6.);
    TH1D h_Mmis_AngleFixcorr2("h_Mmis_AngleFixcorr2", "", 200, -1., 6.);
    TH1D h_Mmis_AngleFixcorr3("h_Mmis_AngleFixcorr3", "", 200, -1., 6.);

    TH2D h_N_mup_mum1("h_N_mup_mum1", "", 11, -0.5, 10.5, 11, -0.5, 10.5);

    TH1D h_Delta_tM1("h_Delta_tM1", "", 200, -1.5, 2);
    TH1D h_Delta_tM_Constrian1("h_Delta_tM_Constrian1", "", 200, -1.5, 2);
    TH1D h_Delta_tM_FixMom1("h_Delta_tM_FixMom1", "", 200, -1.5, 2);
    TH1D h_Delta_tM_FixMom_Constrain1("h_Delta_tM_FixMom_Constrain1", "", 200, -1.5, 2);
    TH1D h_Delta_tM_FixEnergy1("h_Delta_tM_FixEnergy1", "", 200, -1.5, 2);
    TH1D h_Delta_tM_FixEnergy_Constrain1("h_Delta_tM_FixEnergy_Constrain1", "", 200, -1.5, 2);
    TH1D h_Delta_tM_AngleFix1("h_Delta_tM_AngleFix1", "", 200, -1.5, 2);

    TH2D h_Qp2_vs_Q2_Rec1("h_Qp2_vs_Q2_Rec1", "", 200, 0., 10., 200, 0., 16);
    TH2D h_Qp2_vs_Q2_Rec2("h_Qp2_vs_Q2_Rec2", "", 200, 0., 10., 200, 0., 16);
    TH2D h_Qp2_vs_Q2_Rec3("h_Qp2_vs_Q2_Rec3", "", 200, 0., 10., 200, 0., 16);

    TH2D h_Q2_xB_Rec1("h_Q2_xB_Rec1", "", 200, 0., 0.6, 200, 0., 10.);
    TH2D h_Q2_xB_Rec2("h_Q2_xB_Rec2", "", 200, 0., 0.6, 200, 0., 10.);

    TH2D h_xB_tM_Rec1("h_xB_tM_Rec1", "", 200, 0., 0.15 * Eb, 200, 0., 0.75);

    TH2D h_xi_xxGPD1("h_xi_xxGPD1", "", 200, -1., 1., 200, 0., 1.);
    TH2D h_xi_xxGPD2("h_xi_xxGPD2", "", 200, -1., 1., 200, 0., 1.);
    TH2D h_xi_xxGPD7("h_xi_xxGPD7", "", 200, -1., 1., 200, 0., 1.);

    TH2D h_xi_xxGPD_bin_[n_xi_x_bins];
    TH2D h_Qp2_Q2_bin_[n_xi_x_bins];

    for (int i = 0; i < n_xi_x_bins; i++) {
        h_xi_xxGPD_bin_[i] = TH2D(Form("h_xi_xxGPD_bin_%d", i), "", 200, -1., 1., 200, 0., 1.);
        h_Qp2_Q2_bin_[i] = TH2D(Form("h_Qp2_Q2_bin_%d", i), "", 200, 0., 10., 200, 0., 10.);
    }


    TH2D h_xi_xxGPD_MC1("h_xi_xxGPD_MC1", "", 200, -1., 1., 200, 0., 1.);
    TH2D h_xi_xxGPD_MC2("h_xi_xxGPD_MC2", "", 200, -1., 1., 200, 0., 1.);
    TH2D h_xi_xxGPD_MC3("h_xi_xxGPD_MC3", "", 200, -1., 1., 200, 0., 1.);

    TH2D h_dP_P_MC_prot1("h_dP_P_MC_prot1", "", 200, 0., 0.25 * Eb, 200, -0.65, 0.65);
    TH2D h_dE_E_MC_prot1("h_dE_E_MC_prot1", "", 200, 0., 0.25 * Eb, 200, -0.65, 0.65);

    TH1D h_Mmis_MC1("h_Mmis_MC1", "", 200, -1., 6.);
    TH1D h_Mmis_MC2("h_Mmis_MC2", "", 200, -1., 6.);

    TH1D h_Mumu_Corr1("h_Mumu_Corr1", "", 200, 0., 4.);
    TH1D h_Mumu_Corr_AngleFix1("h_Mumu_Corr_AngleFix1", "", 200, 0., 0.35 * Eb);
    TH1D h_Mumu_Corr_AngleFix2("h_Mumu_Corr_AngleFix2", "", 200, 0., 0.35 * Eb);

    TH2D h_xi_xxGPD_Rec1("h_xi_xxGPD_Rec1", "", 200, -0.5, 0.5, 200, 0., 0.5);
    TH2D h_Qp2_Q2_Rec1("h_Qp2_Q2_Rec1", "", 200, 0., 10., 200, 0., 16.);

    TH2D h_MC_th_P_em_xi_x_[n_xi_x_bins][n_xi_x_Q2bins];
    TH2D h_xi_xxGPD_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH1D h_Phi_LH_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D h_Qp2_Q2_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 
    TH2D h_tM_xB_xi_x_[n_xi_x_bins][n_xi_x_Q2bins]; // 

    TH1D h_M_mumu_[n_xi_x_bins];

    for (int i = 0; i < n_xi_x_bins; i++) {

        h_M_mumu_[i] = TH1D(Form("h_M_mumu_%d", i), "", 200, 0., 4.5);
        for (int j = 0; j < n_xi_x_Q2bins; j++) {
            h_xi_xxGPD_xi_x_[i][j] = TH2D(Form("h_xi_xxGPD_xi_x_%d_%d", i, j), "", 200, -0.5, 0.5, 200, 0., 0.5);
            h_Phi_LH_xi_x_[i][j] = TH1D(Form("h_Phi_LH_xi_x_%d_%d", i, j), "", 12, 0, 360);
            h_Qp2_Q2_xi_x_[i][j] = TH2D(Form("h_Qp2_Q2_xi_x_%d_%d", i, j), "", 200, 0., 0.6 * Eb, 200, 0., 0.6 * Eb);
            h_tM_xB_xi_x_[i][j] = TH2D(Form("h_tM_xB_xi_x_%d_%d", i, j), "", 200, 0., 1.2, 200, 0., 1.);
            h_MC_th_P_em_xi_x_[i][j] = TH2D(Form("h_MC_th_P_em_xi_x_%d_%d", i, j), "",  200, 0., 1.05 * Eb, 200, 0., 45.);
        }
    }


    hipo::reader reader;
    //reader.open("Data/Skim_ZeroSuppr_2247_All.hipo");
    reader.open(Form("Data/DDVCS_Run%d_All.hipo", run));

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

            //                        if (evCounter > 200000) {
            //                            break;
            //                        }
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
            int n_em = 0;
            int n_charged = 0;

            vector<RecParticle> v_recp_mum;
            vector<RecParticle> v_recp_mup;
            vector<RecParticle> v_recp_em;
            vector<RecParticle> v_recp_p;

            for (int ipart = 0; ipart < nPart; ipart++) {
                RecParticle recp(bRecPart, bRecCalo, bRecCC, ipart, ind_PCal[ipart], ind_ECin[ipart], ind_ECout[ipart], ind_HTCC[ipart]);

                int nStrip_PCal = recp.duPCal() + recp.dvPCal() + recp.dwPCal();
                int nStrip_ECin = recp.duECin() + recp.dvECin() + recp.dwECin();
                int nStrip_ECout = recp.duECout() + recp.dvECout() + recp.dwECout();

                bool hitAllLayers = (nStrip_PCal >= 3) && (nStrip_ECin >= 3) && (nStrip_ECout >= 3);

                if ((recp.pid() == 211 || recp.pid() == -13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000 && recp.th() < th_MuMax) {

                    bool MIP_Signature = nStrip_PCal <= nMaxPCalStripsPos && nStrip_ECin < nMaxECinStripsPos && nStrip_ECout < nMaxECoutStripsPos;

                    if (!MIP_Signature && hitAllLayers) {
                        continue;
                    }

                    v_recp_mup.push_back(recp);
                    h_th_P_mup1.Fill(recp.p(), recp.th());
                    n_mup = n_mup + 1;
                } else if ((recp.pid() == -211 || recp.pid() == 13) && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000 && recp.th() < th_MuMax) {

                    bool MIP_Signature = nStrip_PCal <= nMaxPCalStripsNeg && nStrip_ECin < nMaxECinStripsNeg && nStrip_ECout < nMaxECoutStripsNeg;

                    if (!MIP_Signature && hitAllLayers) {
                        continue;
                    }

                    v_recp_mum.push_back(recp);
                    h_th_P_mum1.Fill(recp.p(), recp.th());
                    n_mum = n_mum + 1;
                }
            }

            vector<MCPart> v_mc_em; // for merged LUND files we can have more than 1 electron in the Final State
            MCPart mc_em, mc_mum, mc_mup, mc_p;
            for (int imc = 0; imc < nMCPart; imc++) {
                MCPart curPart;

                /*
                 *  For some reasin the MC::Particle bank has several more entries w/ pid == 11 and p = 11 GeV
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

                if (curPart.pid == 11) {

                    mc_em = curPart;
                    v_mc_em.push_back(curPart);

                } else if (curPart.pid == MC_mumPID) {
                    mc_mum = curPart;
                } else if (curPart.pid == MC_mupPID) {
                    mc_mup = curPart;
                    //h_th_P_mupTEST.Fill(L_mup_mc.P(), L_mup_mc.Theta() * TMath::RadToDeg());
                    double p = sqrt(curPart.px * curPart.px + curPart.py * curPart.py + curPart.pz * curPart.pz);
                    double theta = acos(curPart.pz / p);
                    //h_th_P_mupTEST.Fill(p, theta* TMath::RadToDeg());
                } else if (curPart.pid == 2212) {
                    mc_p = curPart;
                }
            }


            /* 
             * Let's find out which e- is inside the acceptance, so we will use that one for further analysis
             */

            for (auto curPart : v_mc_em) {
                TLorentzVector L_em_mc_tmp;
                L_em_mc_tmp.SetPxPyPzE(curPart.px, curPart.py, curPart.pz, curPart.p);
                
                if( emAcc(L_em_mc_tmp) ){
                    mc_em = curPart;
                }
            }

            h_N_mup_mum1.Fill(n_mum, n_mup);

            L_em_mc.SetPxPyPzE(mc_em.px, mc_em.py, mc_em.pz, mc_em.p);
            L_mum_mc.SetPxPyPzE(mc_mum.px, mc_mum.py, mc_mum.pz, sqrt(mc_mum.p * mc_mum.p + Mmu * Mmu));
            L_mup_mc.SetPxPyPzE(mc_mup.px, mc_mup.py, mc_mup.pz, sqrt(mc_mup.p * mc_mup.p + Mmu * Mmu));
            L_p_mc.SetPxPyPzE(mc_p.px, mc_p.py, mc_p.pz, sqrt(mc_p.p * mc_p.p + Mp * Mp));
            h_MC_th_P_em1.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());


            /*
             * Skip events when one of muons goes through the gap between the Moeller cone and the Tungsten shield
             */
            if (L_mum_mc.Theta() * TMath::RadToDeg() < th_MuMin_MC || L_mup_mc.Theta() * TMath::RadToDeg() < th_MuMin_MC) {
                continue;
            }

            TLorentzVector L_q_MC = L_beam - L_em_mc;
            TLorentzVector L_mumu_MC = L_mum_mc + L_mup_mc;
            double Q2_MC = -L_q_MC.M2();
            double Qp2_MC = L_mumu_MC.M2();
            double nue_MC = L_q_MC.E();
            double xB_MC = Q2_MC / (2 * Mp * nue_MC);
            double xiPrime_MC = xB_MC / (2 - xB_MC);
            double xi_MC = xiPrime_MC * (Q2_MC + Qp2_MC) / (Q2_MC);

            //cout<<"Q2 = "<<Q2_MC<<"     Qp2_MC "<<Qp2_MC<<"  nue_MC = "<<nue_MC<<endl;

            TLorentzVector L_mis_MC = L_beam + L_targ - L_mumu_MC - L_em_mc;
            double MMis_mc = L_mis_MC.M2();
            h_Mmis_MC1.Fill(MMis_mc);
            //cout<<MMis_mc<<endl;

            double xx_GPD_MC = 2 * xiPrime_MC - xi_MC; // This is x that GPDs depend on defined as 2*xiPrime  - xi

            h_xi_xxGPD_MC1.Fill(xx_GPD_MC, xi_MC);

            bool em_acc = emAcc(L_em_mc);

            if (em_acc) {
                h_MC_th_P_em_det1.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());
            } else {
                h_MC_th_P_em_miss1.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());
            }
            //cout<<L_em_mc.Px()<<"   "<<L_em_mc.Py()<<
            //cout<<n_mum<<" "<<n_mup<<" "<<em_acc<<endl;

            if (n_mup < 1) {
                h_MC_th_P_mup_miss1.Fill(L_mup_mc.P(), L_mup_mc.Theta() * TMath::RadToDeg());
                h_MC_th_phi_mup_miss1.Fill(L_mup_mc.Phi() * TMath::RadToDeg(), L_mup_mc.Theta() * TMath::RadToDeg());
            } else {
                h_MC_th_P_mup_det1.Fill(L_mup_mc.P(), L_mup_mc.Theta() * TMath::RadToDeg());
                h_MC_th_phi_mup_det1.Fill(L_mup_mc.Phi() * TMath::RadToDeg(), L_mup_mc.Theta() * TMath::RadToDeg());
            }

            if (n_mum < 1) {
                h_MC_th_P_mum_miss1.Fill(L_mum_mc.P(), L_mum_mc.Theta() * TMath::RadToDeg());
                h_MC_th_phi_mum_miss1.Fill(L_mum_mc.Phi() * TMath::RadToDeg(), L_mum_mc.Theta() * TMath::RadToDeg());
            } else {
                h_MC_th_P_mum_det1.Fill(L_mum_mc.P(), L_mum_mc.Theta() * TMath::RadToDeg());
                h_MC_th_phi_mum_det1.Fill(L_mum_mc.Phi() * TMath::RadToDeg(), L_mum_mc.Theta() * TMath::RadToDeg());
            }

            //cout<<L_mum_mc.Theta()<<"   "<<L_mum_mc.Phi()<<"   "<<L_mum_mc.M()<<endl;

            /* 
             * Just mu-mu+ topology
             */
            if( n_mum == 1 && n_mup == 1 ){
                
            }
            
            if (n_mum == 1 && n_mup == 1 && em_acc) {

                RecParticle rec_mup = v_recp_mup.at(0);
                RecParticle rec_mum = v_recp_mum.at(0);

                double th_mup_Rec = rec_mup.th();
                double th_mum_Rec = rec_mum.th();

                double phi_mup_Rec = rec_mup.phi() * TMath::RadToDeg();
                double phi_mum_Rec = rec_mum.phi() * TMath::RadToDeg();

                double th_mup_MC = L_mup_mc.Theta() * TMath::RadToDeg();
                double th_mum_MC = L_mum_mc.Theta() * TMath::RadToDeg();

                double tM_MC = 2 * Mp * (L_p_mc.E() - Mp);

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

                h_MC_th_P_em_FS1.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());
                h_th_P_mup2.Fill(rec_mup.p(), th_mup_Rec);
                h_th_P_mum2.Fill(rec_mum.p(), th_mum_Rec);
                h_MC_th_P_em2.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());

                if ((rec_mup.p() - L_mup_mc.P()) / L_mup_mc.P() > f_Pol3Test2->Eval(rec_mup.p())) {

                    h_th_P_mupTEST1.Fill(rec_mup.p(), th_mup_Rec);
                    h_th_P_mupTEST3.Fill(L_mup_mc.P(), th_mup_MC);
                } else {
                    h_th_P_mupTEST2.Fill(rec_mup.p(), th_mup_Rec);
                    h_th_P_mupTEST4.Fill(L_mup_mc.P(), th_mup_MC);
                }


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

                L_em_rec = L_em_mc;

                L_mup_rec.SetPxPyPzE(rec_mup.px(), rec_mup.py(), rec_mup.pz(), sqrt(rec_mup.p() * rec_mup.p() + Mmu * Mmu));
                L_mum_rec.SetPxPyPzE(rec_mum.px(), rec_mum.py(), rec_mum.pz(), sqrt(rec_mum.p() * rec_mum.p() + Mmu * Mmu));

                TLorentzVector L_mis = L_beam + L_targ - L_em_rec - L_mup_rec - L_mum_rec;

                double Mx2 = L_mis.M2();
                h_Mmis1.Fill(Mx2);

                /*
                 * Let's make kinematic corrections below.
                 */

                double p_mup_corr = L_mup_rec.P() / (1 + f_Elos_mup->Eval(L_mup_rec.P()));
                double momScale = p_mup_corr / L_mup_rec.P();
                L_mup_cor.SetPxPyPzE(L_mup_rec.Px() * momScale, L_mup_rec.Py() * momScale, L_mup_rec.Pz() * momScale, sqrt(p_mup_corr * p_mup_corr + Mmu * Mmu));


                // LorentzVector that has only Momentum corrected, but angles are fixed to generated values
                TLorentzVector L_mup_AngleFixCor = L_mup_cor;
                L_mup_AngleFixCor.SetTheta(L_mup_mc.Theta());
                L_mup_AngleFixCor.SetPhi(L_mup_mc.Phi());


                double p_mum_corr = L_mum_rec.P() / (1 + f_Elos_mum->Eval(L_mum_rec.P()));
                momScale = p_mum_corr / L_mum_rec.P();
                L_mum_cor.SetPxPyPzE(L_mum_rec.Px() * momScale, L_mum_rec.Py() * momScale, L_mum_rec.Pz() * momScale, sqrt(p_mum_corr * p_mum_corr + Mmu * Mmu));
                TLorentzVector L_mis_corr1 = L_beam + L_targ - L_em_rec - L_mup_cor - L_mum_cor;

                // LorentzVector that has only Momentum corrected, but angles are fixed to generated values
                TLorentzVector L_mum_AngleFixCor = L_mum_cor;
                L_mum_AngleFixCor.SetTheta(L_mum_mc.Theta());
                L_mum_AngleFixCor.SetPhi(L_mum_mc.Phi());

                double Mx2Corr = L_mis_corr1.M2();
                h_Mmis_corr1.Fill(Mx2Corr);

                double mup_delta_th_corr = f_ThCorr_mup->Eval(th_mup_Rec);
                RotateTheta(L_mup_cor, mup_delta_th_corr);


                double mum_delta_th_corr = f_ThCorr_mum->Eval(th_mum_Rec);
                RotateTheta(L_mum_cor, mum_delta_th_corr);

                TLorentzVector L_mis_corr2 = L_beam + L_targ - L_em_rec - L_mup_cor - L_mum_cor;
                double Mx2Corr2 = L_mis_corr2.M2();
                h_Mmis_corr2.Fill(Mx2Corr2);


                TLorentzVector L_Mis_Corr_AngleFix = L_beam + L_targ - L_em_rec - L_mup_AngleFixCor - L_mum_AngleFixCor;
                double Mx2CorrAngleFix1 = L_Mis_Corr_AngleFix.M2();
                h_Mmis_AngleFixcorr1.Fill(Mx2CorrAngleFix1);
                h_Mmis_AngleFixcorr_WideRange1.Fill(Mx2CorrAngleFix1);

                h_Mmis_MC2.Fill(MMis_mc);

                TLorentzVector L_MuMu_CorrAngleFix = L_mup_AngleFixCor + L_mum_AngleFixCor;
                double Mmumu = L_MuMu_CorrAngleFix.M();

                if (Mmumu > MinvMin1) {
                    h_Mmis_AngleFixcorr2.Fill(Mx2CorrAngleFix1);
                }
                if (Mmumu > MinvMin2) {
                    h_Mmis_AngleFixcorr3.Fill(Mx2CorrAngleFix1);
                }


                h_Mumu_Corr_AngleFix1.Fill(L_MuMu_CorrAngleFix.M());

                double mup_delta_Phi = -f_PhiCorr_mup->Eval(pt_mup_Rec);
                //cout<<"Before Rotation the angle is "<<L_mup_cor.Theta()*TMath::RadToDeg()<<"  "<<L_mup_cor.Phi()*TMath::RadToDeg()<<endl;
                L_mup_cor.RotateZ(mup_delta_Phi * TMath::DegToRad());
                //cout<<"After Rotation the angle is "<<L_mup_cor.Theta()*TMath::RadToDeg()<<"  "<<L_mup_cor.Phi()*TMath::RadToDeg()<<endl;

                double mum_delta_Phi = -f_PhiCorr_mum->Eval(pt_mum_Rec);
                L_mum_cor.RotateZ(mum_delta_Phi * TMath::DegToRad());

                TLorentzVector L_mis_corr3 = L_beam + L_targ - L_em_rec - L_mup_cor - L_mum_cor;
                double Mx2Corr3 = L_mis_corr3.M2();
                h_Mmis_corr3.Fill(Mx2Corr3);

                h_dP_P_MC_prot1.Fill(L_p_mc.P(), (L_mis_corr3.P() - L_p_mc.P()) / L_p_mc.P());
                h_dE_E_MC_prot1.Fill(L_p_mc.E(), (L_mis_corr3.E() - L_p_mc.E()) / L_p_mc.E());

                TLorentzVector L_mis_Prot_FixMomentum = L_mis_corr3;
                TLorentzVector L_mis_Prot_FixEnergy = L_mis_corr3;

                // In the case of L_mis_Prot_FixMomentum, we will correct the energy, assuming the momentum
                // measurement is more precise
                L_mis_Prot_FixMomentum.SetE(sqrt(L_mis_Prot_FixMomentum.P() * L_mis_Prot_FixMomentum.P() + Mp * Mp));

                double corrected_P = sqrt(L_mis_Prot_FixEnergy.E() * L_mis_Prot_FixEnergy.E() - Mp * Mp);
                double correctionScale_P = corrected_P / L_mis_Prot_FixEnergy.P();
                L_mis_Prot_FixEnergy.SetPxPyPzE(L_mis_Prot_FixEnergy.Px() * correctionScale_P, L_mis_Prot_FixEnergy.Py() * correctionScale_P, L_mis_Prot_FixEnergy.Pz() * correctionScale_P,
                        L_mis_Prot_FixEnergy.E());

                /*
                 *  At this point kinematic corrections are applied.
                 *  We will look into phase space coverage, and estimation of uncertainties.
                 */

                TLorentzVector L_q_Rec = L_beam - L_em_rec; // The LorentzVector of Spacelike photon
                TLorentzVector L_mumu_cor = L_mum_cor + L_mup_cor;

                TVector3 v3_beam_Eprime = L_beam.Vect().Cross(L_em_rec.Vect());
                TVector3 v3_q_qprime = L_q_Rec.Vect().Cross(L_MuMu_CorrAngleFix.Vect());

                double Phi_LH_Rec = L_em_rec.Vect().Dot(v3_q_qprime) > 0 ? v3_beam_Eprime.Angle(v3_q_qprime) * r2d : 360. - v3_beam_Eprime.Angle(v3_q_qprime) * r2d;

                double Q2_Rec = -L_q_Rec.M2();
                double nue_Rec = L_q_Rec.E();

                double Qp2_Rec = Mmumu*Mmumu;

                h_Mumu_Corr1.Fill(Mmumu);

                double xB_Rec = Q2_Rec / (2 * Mp * nue_Rec); // Calculating the xB
                double xiPrime_Rec = xB_Rec / (2 - xB_Rec); // Calculating the xiPrime as xB/(2-xB)
                double xi_Rec = xiPrime_Rec * (Q2_Rec + Qp2_Rec) / (Q2_Rec); // Calculating xi as xiPrime*(Q2 + Qp2)/Q2

                double xx_GPD_Rec = 2 * xiPrime_Rec - xi_Rec; // This is x that GPDs depend on defined as 2*xiPrime - xi

                double tM_Rec = -(L_mis_corr3 - L_targ).M2();
                double tM_AngleFix = -(L_Mis_Corr_AngleFix - L_targ).M2();

                double tM_Rec_FixMom = -(L_mis_Prot_FixMomentum - L_targ).M2();
                double tM_Rec_FixMomConstrain = 2 * Mp * (L_mis_Prot_FixMomentum.E() - Mp);

                double tM_Rec_FixEnergy = -(L_mis_Prot_FixEnergy - L_targ).M2();
                double tM_Rec_FixEnergyConstrain = 2 * Mp * (L_mis_Prot_FixEnergy.E() - Mp);

                /*
                 *  We constrain below the mass of the missing particle to the
                 * proton mass. We should expect better resolution with this 
                 * way of calculation of tM
                 */
                double tM_RecConstrain = 2 * Mp * (sqrt(L_mis_corr3.P() * L_mis_corr3.P() + Mp * Mp) - Mp);

                h_Delta_tM_AngleFix1.Fill(tM_AngleFix - tM_MC);
                h_Delta_tM1.Fill(tM_Rec - tM_MC);
                h_Delta_tM_Constrian1.Fill(tM_RecConstrain - tM_MC);

                h_Delta_tM_FixMom1.Fill(tM_Rec_FixMom - tM_MC);
                h_Delta_tM_FixMom_Constrain1.Fill(tM_Rec_FixMomConstrain - tM_MC);
                h_Delta_tM_FixEnergy1.Fill(tM_Rec_FixEnergy - tM_MC);
                h_Delta_tM_FixEnergy_Constrain1.Fill(tM_Rec_FixEnergyConstrain - tM_MC);

                //cout<<tM_Rec_FixEnergyConstrain<< "   " << tM_Rec_FixEnergy << "   "<< tM_Rec<<"   " << tM_Rec_FixEnergyConstrain - tM_Rec_FixEnergy<<endl;

                h_Qp2_vs_Q2_Rec1.Fill(Q2_Rec, Qp2_Rec);

                h_xB_tM_Rec1.Fill(tM_AngleFix, xB_Rec);
                h_Q2_xB_Rec1.Fill(xB_Rec, Q2_Rec);

                h_xi_xxGPD1.Fill(xx_GPD_Rec, xi_Rec);
                h_xi_xxGPD_MC2.Fill(xx_GPD_MC, xi_MC);


                /*
                 * Applying Missing Mass cut
                 */
                if (Mx2CorrAngleFix1 < MM2_min || Mx2CorrAngleFix1 > MM2_max) {
                    continue;
                }

                h_MC_th_P_em_FS2.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());

                //cout<<xx_GPD_MC<<"   "<<xi_MC<<endl;
                h_Q2_xB_Rec2.Fill(xB_Rec, Q2_Rec);
                h_Mumu_Corr_AngleFix2.Fill(L_MuMu_CorrAngleFix.M());

                int x_xi_Bin = findX_Xi_Bin(xx_GPD_Rec, xi_Rec);

                if (x_xi_Bin >= 0) {
                    h_xi_xxGPD_bin_[x_xi_Bin].Fill(xx_GPD_Rec, xi_Rec);
                    h_Qp2_Q2_bin_[x_xi_Bin].Fill(Q2_Rec, Qp2_Rec);
                }

                if (xx_GPD_MC < bin0_x_max && xx_GPD_MC > bin0_x_min && xi_MC > bin0_xi_min && xi_MC < bin0_xi_max) {
                    h_xi_xxGPD2.Fill(xx_GPD_Rec, xi_Rec);
                    h_xi_xxGPD_MC3.Fill(xx_GPD_MC, xi_MC);
                }

                if (tM_AngleFix < t_Max_xi_xiStudy) {
                    // cout << "Kuku" << endl;

                    h_MC_th_P_em_FS3.Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());
                    
                    h_xi_xxGPD_Rec1.Fill(xx_GPD_Rec, xi_Rec);
                    h_Qp2_Q2_Rec1.Fill(Q2_Rec, Qp2_Rec);

                    int bin_xi_x = -1;
                    int bin_Q2 = -1;

                    if (xx_GPD_Rec > bin0_x_min && xx_GPD_Rec < bin0_x_max && xi_Rec > bin0_xi_min && xi_Rec < bin0_xi_max) {
                        bin_xi_x = 0;
                        bin_Q2 = Qp2_Rec > Qp2Q2_bin0_Offset + Q2_Rec * Qp2Q2_bin0_Slope ? 1 : 0;
                        //bin_Q2 = Qp2_Rec > 4.6 - Q2_Rec / 1.65 ? 1 : 0;
                    } else if (xx_GPD_Rec > bin1_x_min && xx_GPD_Rec < bin1_x_max && xi_Rec > bin1_xi_min && xi_Rec < bin1_xi_max) {
                        bin_xi_x = 1;
                        bin_Q2 = Qp2_Rec > Qp2Q2_bin1_Offset + Q2_Rec * Qp2Q2_bin1_Slope ? 1 : 0;
                        //bin_Q2 = Qp2_Rec > 6.9 - Q2_Rec / 0.55 ? 1 : 0;
                    }

                    //cout<<bin_xi_x<<"   "<<bin_Q2<<endl;

                    if (bin_xi_x >= 0 && bin_Q2 >= 0) {
                        h_M_mumu_[bin_xi_x].Fill(L_MuMu_CorrAngleFix.M());
                    }

                    if (bin_xi_x >= 0 && bin_Q2 >= 0) {
                        h_Phi_LH_xi_x_[bin_xi_x][bin_Q2].Fill(Phi_LH_Rec);
                        h_xi_xxGPD_xi_x_[bin_xi_x][bin_Q2].Fill(xx_GPD_Rec, xi_Rec);
                        h_Qp2_Q2_xi_x_[bin_xi_x][bin_Q2].Fill(Q2_Rec, Qp2_Rec);
                        h_tM_xB_xi_x_[bin_xi_x][bin_Q2].Fill(tM_AngleFix, xB_Rec);
                        h_MC_th_P_em_xi_x_[bin_xi_x][bin_Q2].Fill(L_em_mc.P(), L_em_mc.Theta() * TMath::RadToDeg());
                    }


                }


            }

        }
    } catch (const char* msg) {
        cerr << msg << endl;
    }

    //double weight = xsec[0] * pbn * totLumi / nev;
    double weight = tot_xSec * pbn * totLumi / evCounter;

    //double weight = 1.;

    TList *l_objs = gDirectory->GetList();

    for (TIter curObj = l_objs->begin(); curObj != l_objs->end(); curObj.Next()) {
        ((TH1*) * curObj)->Scale(weight);
    }

    cout << "\n\n" << endl;
    cout << "Weight = " << weight << endl;
    gDirectory->Write();
    return 0;
}

bool emAcc(TLorentzVector& L) {

    const double th_max = 29.5; // deg
    const double th_min = 7.5; // deg
    const double P_max = 22.; // GeV
    const double P_min = 1.; // GeV

    return L.Theta() * TMath::RadToDeg() > th_min && L.Theta() * TMath::RadToDeg() < th_max && L.P() > P_min && L.P() < P_max;
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

int findX_Xi_Bin(double x, double xi) {
    const double bin0_x_max = -0.05;
    const double bin0_x_min = -0.09;
    const double bin0_xi_max = 0.245;
    const double bin0_xi_min = 0.225;

    const double bin1_x_max = 0.09;
    const double bin1_x_min = 0.05;
    const double bin1_xi_max = 0.245;
    const double bin1_xi_min = 0.225;

    int bin = -1;

    if (x > bin0_x_min && x < bin0_x_max && xi > bin0_xi_min && xi < bin0_xi_max) {
        bin = 0;
    } else if (x > bin1_x_min && x < bin1_x_max && xi > bin1_xi_min && xi < bin1_xi_max) {
        bin = 1;
    }

    return bin;
}
