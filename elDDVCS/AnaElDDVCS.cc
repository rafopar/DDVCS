//
// Created by rafopar on 5/22/26.
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
#include <RecEvent.h>

#include "DDVCSTools.h"
#include "include/AnaElDDVCSHelper.h"

#include  <AnaDB.h>
#include <ATen/core/interned_strings.h>

#include "ElectronDDVCSKine.h"

using namespace std;

int main(const int argc, char *argv[]) {
    if (argc < 3) {
        cout << "Usage: AnaElDDVCS.exe <run> <dataset>" << endl;
        cout << "  dataset: 8p477GeV or 6p395GeV" << endl;
        exit(1);
    }

    int run = atoi(argv[1]);
    string dataset = argv[2];
    char inputFile[256];

    sprintf(inputFile, "SkimDDVCS/%s/elDDVCS_%06d.hipo", dataset.c_str(), run);

    // Reading the DB
    AnaDB db("/work/clas12/rafopar/DB/RafoAnaDB.db");

    static constexpr double vlight = 29.373;

    // ------- INIT Analysis cuts -------
    AnaCuts cuts;
    InitCuts(db, run, cuts);
    TFile file_out(Form("OutHists/AnaElDDVCS/Hists_%d.root", run), "Recreate");

    HistoList *hists = new HistoList;
    InitHistos(hists);

    ElectronDDVCSKine ddvcs_kine_NoProt;

    ddvcs_kine_NoProt.SetEb(cuts.Eb);
    ddvcs_kine_NoProt.SetMtarg(DDVCSTools::Mprot);



    hipo::reader reader;
    reader.open(inputFile);
    hipo::dictionary dict;
    reader.readDictionary(dict);

    hipo::bank bScint(dict.getSchema("REC::Scintillator"));
    hipo::bank bRunConf(dict.getSchema("RUN::config"));

    RecEvent recEvent(dict);
    hipo::event event;

    try {
        while (reader.next()) {
            reader.read(event);
            recEvent.Load(event);

             std::vector<const RecPart*> v_el, v_ep, v_prot;

            for (auto &p: recEvent.Particles()) {
                // cout<<"PID = "<<p.pid()<<endl;
                //    if (p.PCal().valid)  cout << p.PCal().energy<<endl;
                //    if (p.HTCC().valid)  cout << p.HTCC().nphe<<endl;

                if (p.pid() == 11) {
                    hists->h_chi2PID_em1->Fill(p.chi2pid());
                    hists->h_vz_em1->Fill(p.vz());

                    double vt = p.FTOF1B().time - p.FTOF1B().path/vlight - recEvent.starttime();
                    hists->h_vt_em1_1->Fill(vt);

                    bool isPCalEmin = p.PCal().energy > cuts.PCalEMin;
                    bool isPCalVmin = p.PCal().lv > cuts.PCalVmin;
                    bool isPCalWmin = p.PCal().lw > cuts.PCalWmin;
                    bool isVz = p.vz() > cuts.emVzMin && p.vz() < cuts.emVzMax;

                    bool isVt = vt > cuts.emVtMin && vt < cuts.emVtMax;

                    if ( isPCalEmin && isPCalVmin && isPCalWmin && isVz && isVt ) {
                        v_el.push_back(&p);
                    }


                }else if (p.pid() == -11) {

                    double Etot = p.PCal().energy + p.ECin().energy + p.ECin().energy;

                    hists->h_vz_ep1->Fill(p.vz());
                    hists->h_chi2PID_ep1->Fill(p.chi2pid());
                    hists->h_EPCal_ep1->Fill(p.PCal().energy);
                    hists->h_EECin_PCal_ep1->Fill(p.PCal().energy, p.ECin().energy);
                    hists->h_SF_PCal_ep1->Fill(p.p(), p.PCal().energy/p.p());
                    hists->h_SF_Tot_ep1->Fill(p.p(),  Etot/p.p());
                    hists->h_SF_ECin_PCal_ep1->Fill(p.PCal().energy/p.p(), p.ECin().energy/p.p());

                    double vt = p.FTOF1B().time - p.FTOF1B().path/vlight - recEvent.starttime();
                    hists->h_vt_ep1->Fill(vt);

                    bool isPCalEmin = p.PCal().energy > cuts.PCalEMin;
                    bool isPCalVmin = p.PCal().lv > cuts.PCalVmin;
                    bool isPCalWmin = p.PCal().lw > cuts.PCalWmin;
                    bool isvz = p.vz() > cuts.epVzMin && p.vz() < cuts.epVzMax;

                    bool isVt = vt > cuts.epVtMin && vt < cuts.epVtMax;

                    if ( isPCalEmin && isPCalVmin && isPCalWmin && isvz && isVt ) {
                        v_ep.push_back(&p);
                    }
                }else if ( p.pid() == 2212 ) {
                    hists->h_vz_prot1->Fill(p.vz());
                    hists->h_chi2PID_prot1->Fill(p.chi2pid());

                    if ( p.chi2pid() > cuts.protChi2PIDMin && p.chi2pid() < cuts.protChi2PIDMax ) {
                        if ( !p.FTOF1B().valid && !p.FTOF1A().valid && !p.FTOF2().valid && !p.CTOF().valid ) {
                            cout<<"No any scinitillator hit..."<<endl;
                        }

                        const ScintillatorResponse *SC_prot;
                        if (p.FTOF1B().valid) {
                            SC_prot = &p.FTOF1B();
                        }else if ( p.FTOF1A().valid ) {
                            SC_prot = &p.FTOF1A();
                        }else if (p.FTOF2().valid ) {
                            SC_prot = &p.FTOF2();
                        }else if (p.CTOF().valid ) {
                            SC_prot = &p.CTOF();
                        }

                        //double vt = SC_prot->time - SC_prot->path/(p.beta()*vlight) - recEvent.starttime();
                        double vt = SC_prot->time - SC_prot->path/((p.p()/sqrt(p.p()*p.p() + DDVCSTools::Mprot*DDVCSTools::Mprot ))*vlight) - recEvent.starttime();
                        hists->h_vt_prot1->Fill(vt);
                        v_prot.push_back(&p);
                    }
                }
            }

            if ( v_el.size() == 2 && v_ep.size() == 1 ) {
                ROOT::Math::PxPyPzEVector L_em1(v_el.at(0)->px(), v_el.at(0)->py(), v_el.at(0)->pz(), v_el.at(0)->p());
                ROOT::Math::PxPyPzEVector L_em2(v_el.at(1)->px(), v_el.at(1)->py(), v_el.at(1)->pz(), v_el.at(1)->p());
                ROOT::Math::PxPyPzEVector L_ep(v_ep.at(0)->px(), v_ep.at(0)->py(), v_ep.at(0)->pz(), v_ep.at(0)->p());

                int n_prot = v_prot.size();

                // double vt_el1 = v_el.at(0)->FTOF1B().time - v_el.at(0)->FTOF1B().path/vlight - recEvent.starttime();
                // double vt_el2 = v_el.at(1)->FTOF1B().time - v_el.at(1)->FTOF1B().path/vlight - recEvent.starttime();
                //
                // hists->h_vt_em1_1->Fill(vt_el1);
                // hists->h_vt_em2_1->Fill(vt_el2);


                ddvcs_kine_NoProt.SetKineem1em2ep(L_em1, L_em2, L_ep);

                double MxRecoil = ddvcs_kine_NoProt.GetMxRecoil();

                hists->h_MXRecoil1->Fill(MxRecoil);
                if ( n_prot == 1 ) {
                    hists->h_MXRecoil2->Fill(MxRecoil);
                }

                double Minv1 = ddvcs_kine_NoProt.GetMinv_1();
                double Minv2 = ddvcs_kine_NoProt.GetMinv_2();

                hists->h_Minv12_1->Fill(Minv1, Minv2);
                hists->h_MxRecoil_Minv1_1->Fill(Minv1, MxRecoil);
                hists->h_MxRecoil_Minv2_1->Fill(Minv2, MxRecoil);

                if ( Minv1 > 1.1 || Minv2 > 1.1 ) {
                    hists->h_MXRecoil3->Fill(MxRecoil);
                }
            }
        }
    } catch (std::exception &e) {
        cerr << e.what() << endl;
    }

    gDirectory->Write();
    file_out.Close();
}
