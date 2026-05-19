//
// Created by rafopar on 4/8/26.
//

#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
#include <RecParticle.h>

#include "DDVCSTools.h"
#include "include/elDDVCSHelper.h"

#include  <AnaDB.h>
#include <ATen/core/interned_strings.h>

#include "ElectronDDVCSKine.h"


using namespace std;

int main(const int argc, char* argv[]) {

        if (argc < 2) {
        cout << "Please provide the run number." << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    int run = atoi(argv[1]);
    char inputFile[256];

    sprintf(inputFile, "Train/skim4_%06d.hipo", run);

    TFile file_out(Form("SkimDDVCS/SkimDDVCS_%d.root", run), "Recreate");

    DDVCSTools::Eb = 8.477; ////// FIX ME, THIS SHOULD BE READOUT FROM THE DB //////

    ElectronDDVCSKine ddvcs_kine_NoProt;
    ddvcs_kine_NoProt.SetEb(DDVCSTools::Eb);
    ddvcs_kine_NoProt.SetMtarg(DDVCSTools::Mprot);


    // Reading the DB
    AnaDB db("/work/clas12/rafopar/DB/RafoAnaDB.db");

    HistoList *hists = new HistoList;
    InitHistos(hists);


    // ------- INIT Analysis cuts -------
    AnaCuts cuts;
    InitCuts(db, run, cuts);

    cout<<cuts.emVzMin<<endl;

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

    constexpr int nMax_samePID{20};
    int ind_em[nMax_samePID];
    int ind_ep[nMax_samePID];
    int ind_prot[nMax_samePID];
    int sec_em[nMax_samePID];
    int sec_ep[nMax_samePID];

    int evCounter = 0;
    try {
        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            // if (evCounter > 40000) {
            //     break;
            // }

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
                    if (bRecCC.getInt("detector", i_cc) == static_cast<int>(DDVCSTools::detType::DET_HTCC)) {
                        if (bRecCC.getInt("pindex", i_cc) == i_part) {
                            ind_HTCC[i_part] = i_cc;
                        }
                    }
                }

                int nTOF = bRecSC.getRows();

                // ===================== FTOF ========================
                for (int itof = 0; itof < nTOF; itof++) {
                    if (bRecSC.getInt("detector", itof) == static_cast<int>(DDVCSTools::detType::DET_FTOF) && bRecSC.getInt("layer", itof) == static_cast<int>(DDVCSTools::detType::panel1b)) {
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

                        if (layer == static_cast<int>(DDVCSTools::detType::PCal_Layer)) {
                            ind_PCal[i_part] = i_cal;
                        } else if (layer == static_cast<int>(DDVCSTools::detType::ECin_Layer)) {
                            ind_ECin[i_part] = i_cal;
                        } else if (layer == static_cast<int>(DDVCSTools::detType::ECout_Layer)) {
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

                if (recp.charge() != 0) {
                    n_charged++;
                }

                double Etot = recp.energyPCal() + recp.energyECin() + recp.energyECout();
                if (recp.pid() == 11 && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000) {
                    hists->h_vz_em1->Fill(recp.vz());
                    hists->h_chi2PID_em1->Fill(recp.chi2pid());
                    hists->h_EPCal_em1->Fill(recp.energyPCal());
                    hists->h_EECin_PCal_em1->Fill(recp.energyPCal(), recp.energyECin());
                    hists->h_SF_PCal_em1->Fill(recp.p(), recp.energyPCal()/recp.p());
                    hists->h_SF_Tot_em1->Fill(recp.p(), Etot/recp.p());
                    hists->h_SF_ECin_PCal_em1->Fill(recp.energyPCal()/recp.p(), recp.energyECin()/recp.p());

                    bool isPCalEmin = recp.energyPCal() > cuts.PCalEMin;
                    bool isPCalVmin = recp.lvPCal() > cuts.PCalVmin;
                    bool isPCalWmin = recp.lvPCal() > cuts.PCalWmin;
                    bool isVz = recp.vz() > cuts.emVzMin && recp.vz() < cuts.emVzMax;

                    if ( isPCalEmin && isPCalVmin && isPCalWmin && isVz ) {
                        ind_em[n_em] = ipart;
                        sec_em[n_em] = recp.phi() / 60;
                        n_em = n_em + 1;
                    }


                }else if ( recp.pid() == -11 && TMath::Abs(recp.status()) >= 2000 && TMath::Abs(recp.status()) < 4000 ) {
                    hists->h_vz_ep1->Fill(recp.vz());
                    hists->h_chi2PID_ep1->Fill(recp.chi2pid());
                    hists->h_EPCal_ep1->Fill(recp.energyPCal());
                    hists->h_EECin_PCal_ep1->Fill(recp.energyPCal(), recp.energyECin());
                    hists->h_SF_PCal_ep1->Fill(recp.p(), recp.energyPCal()/recp.p());
                    hists->h_SF_Tot_ep1->Fill(recp.p(), Etot/recp.p());
                    hists->h_SF_ECin_PCal_ep1->Fill(recp.energyPCal()/recp.p(), recp.energyECin()/recp.p());

                    bool isPCalEmin = recp.energyPCal() > cuts.PCalEMin;
                    bool isPCalVmin = recp.lvPCal() > cuts.PCalVmin;
                    bool isPCalWmin = recp.lvPCal() > cuts.PCalWmin;
                    bool isvz = recp.vz() > cuts.epVzMin && recp.vz() < cuts.epVzMax;

                    if ( ep_PCalECin_SF_Cut(cuts, recp) && isPCalEmin && isPCalVmin && isPCalWmin && isvz ) {
                        ind_ep[n_ep] = ipart;
                        sec_ep[n_ep] = recp.phi() / 60;
                        n_ep = n_ep + 1;
                    }

                }else if (recp.pid() == 2212) {
                    hists->h_vz_prot1->Fill(recp.vz());
                    hists->h_chi2PID_prot1->Fill(recp.chi2pid());

                    if ( recp.chi2pid() > cuts.protChi2PIDMin && recp.chi2pid() < cuts.protChi2PIDMax ) {
                        ind_prot[n_prot] = ipart;
                        n_prot = n_prot + 1;
                    }

                }

            }

            hists->h_n_ep_em1->Fill(n_em, n_ep);
            hists->h_n_em1->Fill(n_em);
            hists->h_n_ep1->Fill(n_ep);
            hists->h_n_p1->Fill(n_prot);


            if ( n_em == 2 && n_ep == 1 ) {

                TLorentzVector L_em1, L_em2, L_ep;

                RecParticle part_em1(bRecPart, bRecCalo, bRecCC, ind_em[0], ind_PCal[ind_em[0]], ind_ECin[ind_em[0]],
                     ind_ECout[ind_em[0]], ind_HTCC[ind_em[0]]);

                RecParticle part_em2(bRecPart, bRecCalo, bRecCC, ind_em[1], ind_PCal[ind_em[1]], ind_ECin[ind_em[1]],
                                     ind_ECout[ind_em[1]], ind_HTCC[ind_em[1]]);
                RecParticle part_ep(bRecPart, bRecCalo, bRecCC, ind_ep[0], ind_PCal[ind_ep[0]], ind_ECin[ind_ep[0]],
                                    ind_ECout[ind_ep[0]], ind_HTCC[ind_ep[0]]);

                L_em1.SetPxPyPzE(part_em1.px(), part_em1.py(), part_em1.pz(), part_em1.p());
                L_em2.SetPxPyPzE(part_em2.px(), part_em2.py(), part_em2.pz(), part_em2.p());
                L_ep.SetPxPyPzE(part_ep.px(), part_ep.py(), part_ep.pz(), part_ep.p());
                ddvcs_kine_NoProt.SetKineem1em2ep(&L_em1, &L_em2, &L_ep);

                double MxRecoil = ddvcs_kine_NoProt.GetMxRecoil();

                hists->h_MXRecoil1->Fill(MxRecoil);

                double m_emep1 = ddvcs_kine_NoProt.GetMinv_1();
                double m_emep2 = ddvcs_kine_NoProt.GetMinv_2();

                hists->h_Minv12_1->Fill(m_emep1, m_emep2);

                if ( m_emep1 > 0.25 && m_emep2 > 0.25 ) {
                    hists->h_MXRecoil2->Fill(MxRecoil);
                }

            }


        }
    }catch (const char *msg) {
        cerr << msg << endl;
    }

    gDirectory->Write();
    file_out.Close();
    return 0;
}