//
// Created by rafopar on 12/18/25.
//

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <iostream>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
#include <TLorentzVector.h>
using namespace std;

bool isMIP( hipo::bank &bRecCalo, int ind_PCal, int ind_ECin, int ind_ECout);

int main( int argc, char* argv[] ) {

    if ( argc != 2 ) {
        cerr << "Usage: DebugSkim <Run>" << endl;
        cerr << "Exiting..." << endl;
        exit( 1 );
    }

    int run = atoi( argv[1] );

    char inputFile[256];

    const double MOM_HIGH = 2.;
    const double Mp = 0.9383;
    const double Mmu = 0.1057;
    const double Mpip = 0.13957;
    const double Mpim = 0.13957;

    auto file_out = new TFile(Form("DebugSkim_%d.root", run), "Recreate");

    TH2D h_P_pim_pip("h_P_pim_pip", "", 200, 0., 5., 200, 0., 5.);
    TH2D h_P_pim_prot("h_P_pim_prot", "", 200, 0., 5., 200, 0., 5.);
    TH2D h_P_mum_mup("h_P_mum_mup", "", 200, 0., 5., 200, 0., 5.);
    TH2D h_P_mum_prot("h_P_mum_prot", "", 200, 0., 5., 200, 0., 5.);
    TH1D h_Delta_t_pip_P1("h_Delta_t_pip_P1", "", 200, -1., 5.);
    TH2D h_Delta_t_vs_delta_phi_P1("h_Delta_t_vs_delta_phi_P1", "", 200, -1., 5., 200, -180., 180.);

    sprintf(inputFile, "Data/TCSJPsiSkim/jpsitcs_00%d.hipo", run);

    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;

    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));


    int evCounter = 0;

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

                std::vector<TLorentzVector> v_L_electrons;
                std::vector<TLorentzVector> v_L_FTelectrons;
                std::vector<TLorentzVector> v_L_positrons;
                std::vector<TLorentzVector> v_L_protons;

                std::vector<TLorentzVector> v_L_mups;
                std::vector<TLorentzVector> v_L_mums;
                std::vector<TLorentzVector> v_L_pips;
                std::vector<TLorentzVector> v_L_pims;

                map<int, int> ind_PCal;
                map<int, int> ind_ECin;
                map<int, int> ind_ECout;


                int nPart = bRecPart.getRows();

                for (int i_part = 0; i_part < nPart; i_part++) {
                    ind_PCal[i_part] = -1;
                    ind_ECin[i_part] = -1;
                    ind_ECout[i_part] = -1;

                    int nCal = bRecCalo.getRows();
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

                int npos = 0;
                int nposFD = 0;
                int nneg = 0;
                int ncharge = 0;

                int nFT_any = 0;

                for (int ip = 0; ip < nPart; ip++) {

                    int pid = bRecPart.getInt("pid", ip);

                    double px = bRecPart.getFloat("px", ip);
                    double py = bRecPart.getFloat("py", ip);
                    double pz = bRecPart.getFloat("pz", ip);
                    int charge = bRecPart.getInt("charge", ip);
                    int status = bRecPart.getInt("status", ip);

                    bool isFT = TMath::Abs(status) < 2000;
                    bool isFD = TMath::Abs(status) >= 2000 && TMath::Abs(status) < 4000;

                    if ( isFT ) {
                        nFT_any = nFT_any + 1;
                    }

                    if ( charge == 0 ){continue;}


                    if (charge > 0) {
                        npos = npos + 1;
                        if ( isFD ) {
                            nposFD = nposFD + 1;
                        }
                    }else {
                        nneg = nneg + 1;
                    }
                    ncharge = ncharge + 1;


                    if (pid == 11 && isFD) {
                        v_L_electrons.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz));
                    }else if (pid == -11 && isFD) {
                        v_L_positrons.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz));
                    }else if (pid == 2212 && isFD) {
                        v_L_protons.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz + Mp*Mp));
                    }else if (pid == 11 && isFT) {
                        v_L_FTelectrons.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz));
                    }else if (pid == 211 && isFD) {
                        v_L_pips.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz + Mpip*Mpip));
                    }else if (pid == -211 && isFD ) {
                        v_L_pims.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz + Mpim*Mpim));
                    }



                    if ( isMIP(bRecCalo, ind_PCal[ip], ind_ECin[ip], ind_ECout[ip]) ) {
                        if (charge > 0) {
                            v_L_mups.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz + Mmu*Mmu));
                        }else {
                            v_L_mums.emplace_back(px, py, pz, sqrt(px * px + py * py + pz * pz + Mmu*Mmu));
                        }

                    }
                }

                // Check Valery's Exclusive condition: which is 1pi+, 1pi- and 1 proton

                if ( v_L_pips.size() == 1 && v_L_pims.size() == 1 && v_L_protons.size() == 1 && ncharge == 3 && nFT_any  == 0) {
                    h_P_pim_pip.Fill(v_L_pips.at(0).P(), v_L_pims.at(0).P() );
                    h_P_pim_prot.Fill(v_L_protons.at(0).P(), v_L_pims.at(0).P());

                    h_Delta_t_pip_P1.Fill( v_L_pips.at(0).P() - v_L_protons.at(0).P() );

                    h_Delta_t_vs_delta_phi_P1.Fill(  v_L_pips.at(0).P() - v_L_protons.at(0).P(), (v_L_pips.at(0).Phi() - v_L_protons.at(0).Phi())*TMath::RadToDeg() );

                    // if ( v_L_pips.at(0).P() - v_L_protons.at(0).P() < 0 ) {
                    //     bRecPart.show();
                    // }


                    // if ( v_L_pims.at(0).P() < 0.8 ) {
                    //     cout<<"n_mum = "<<v_L_mums.size()<<"  n_mup = "<<v_L_mums.size()<<"   n_pos = "<<npos<<"  nposFD = "<<nposFD<<"   n_em = "<<v_L_electrons.size()<<"   n_ep = "<<v_L_positrons.size()<<"   n_em_FT = "<<v_L_FTelectrons.size()<<"  nneg = "<<nneg<<endl;
                    //     cout<<"indPCal 0 = "<<ind_PCal[2]<<"   indECin 0 = "<<ind_ECin[2]<<"    ind_ECout 0 = "<<ind_ECout[2]<<endl;
                    //     bRecPart.show();
                    //     bRecCalo.show();
                    // }
                }
                if ( v_L_mums.size() == 1 && v_L_mups.size() == 1 && v_L_protons.size() == 1 && ncharge == 3 ) {
                    h_P_mum_mup.Fill( v_L_mups.at(0).P(), v_L_mums.at(0).P() );
                    h_P_mum_prot.Fill(v_L_protons.at(0).P(), v_L_mums.at(0).P());
                }

            }
        }catch (exception& e) {
            cerr << "Exception: " << e.what() << endl;
        }

    gDirectory->Write();
    return 0;
}

bool isMIP( hipo::bank &bRecCalo, int ind_PCal, int ind_ECin, int ind_ECout) {
    const double MIP_ECOUT_MAX = 0.11;
    const double MIP_ECIN_MAX = 0.1;
    const double MIP_PCAL_MAX = 0.2;

    if ( ind_PCal < 0 || ind_ECin < 0 || ind_ECout < 0 ) {return false;}

    double E_PCal = bRecCalo.getFloat("energy", ind_PCal);
    double E_ECin = bRecCalo.getFloat("energy", ind_ECin);
    double E_ECOut = bRecCalo.getFloat("energy", ind_ECout);

    if ( E_PCal < 0 || E_ECin < 0 || E_ECOut < 0 ) {
        return false;
    }

    return E_PCal < MIP_PCAL_MAX && E_ECin < MIP_ECIN_MAX && E_ECOut < MIP_ECOUT_MAX;
}
