//
// Created by rafopar on 5/25/26.
//

#ifndef DDVCS_STUDIES_ANAELDDVCSHELPER_H
#define DDVCS_STUDIES_ANAELDDVCSHELPER_H

#include <TH1.h>
#include <AnaDB.h>
//#include <VecCore/Backend/Implementation.h>

struct AnaCuts{
    double Eb;

    double emVzMin;
    double emVzMax;
    double emVtMin;
    double emVtMax;
    double epVzMin;
    double emPCalSFMin;
    double emSFSumMin;
    double epVzMax;
    double epVtMin;
    double epVtMax;
    double epPCalSFMin;
    double epSFSumMin;
    double protVzMin;
    double protVzMax;
    double protVtMin;
    double protVtMax;

    double PCalEMin;
    double PCalVmin;
    double PCalWmin;
    double PCalUmin;

    double emChi2PIDMin;
    double emChi2PIDMax;
    double epChi2PIDMin;
    double epChi2PIDMax;
    double protChi2PIDMin;
    double protChi2PIDMax;
};

struct HistoList {

    /// ------------------------- PID histograms -------------------------------------
    TH1D *h_vz_em1;
    TH1D *h_vz_ep1;
    TH1D *h_vz_prot1;
    TH1D *h_chi2PID_em1;
    TH1D *h_chi2PID_ep1;
    TH1D *h_chi2PID_prot1;

    TH1D *h_EPCal_em1;
    TH1D *h_EPCal_ep1;

    TH2D *h_EECin_PCal_em1;
    TH2D *h_EECin_PCal_ep1;

    TH2D *h_SF_PCal_em1;
    TH2D *h_SF_PCal_ep1;

    TH2D *h_SF_Tot_em1;
    TH2D *h_SF_Tot_ep1;

    TH2D *h_SF_ECin_PCal_em1;
    TH2D *h_SF_ECin_PCal_ep1;


    /// ------------------------- FS histograms -------------------------------------
    TH1D *h_MXRecoil1;
    TH1D *h_MXRecoil2;
    TH1D *h_MXRecoil3;

    TH1D *h_vt_em1_1;
    TH1D *h_vt_em2_1;
    TH1D *h_vt_ep1;
    TH1D *h_vt_prot1;

    TH2D *h_Minv12_1;
    TH2D *h_MxRecoil_Minv1_1;
    TH2D *h_MxRecoil_Minv2_1;
};

void InitHistos(HistoList *hists) {
    hists->h_vz_em1 = new TH1D("h_vz_em1", "", 200, -25., 35.);
    hists->h_vz_ep1 = new TH1D("h_vz_ep1", "", 200, -25., 35.);
    hists->h_vz_prot1 = new TH1D("h_vz_prot1", "", 200, -25., 35.);
    hists->h_chi2PID_em1 = new TH1D("h_chi2PID_em1", "", 200, -8., 8.);
    hists->h_chi2PID_ep1 = new TH1D("h_chi2PID_ep1", "", 200, -8., 8.);
    hists->h_chi2PID_prot1 = new TH1D("h_chi2PID_prot1", "", 200, -50., 50.);

    hists->h_EPCal_em1 = new TH1D("h_EPCal_em1", "", 200, 0., 2.);
    hists->h_EPCal_ep1 = new TH1D("h_EPCal_ep1", "", 200, 0., 2.);

    hists->h_EECin_PCal_em1 = new TH2D("h_EECin_PCal_em1", "", 200, 0., 2., 200, 0., 2.);
    hists->h_EECin_PCal_ep1 = new TH2D("h_EECin_PCal_ep1", "", 200, 0., 2., 200, 0., 2.);

    hists->h_SF_PCal_em1 = new TH2D("h_SF_PCal_em1", "", 200, 0., 0.95*DDVCSTools::Eb, 200, 0., 0.5);
    hists->h_SF_PCal_ep1 = new TH2D("h_SF_PCal_ep1", "", 200, 0., 0.95*DDVCSTools::Eb, 200, 0., 0.5);

    hists->h_SF_Tot_em1 = new TH2D("h_SF_Tot_em1", "", 200, 0., 0.95*DDVCSTools::Eb, 200, 0., 0.5);
    hists->h_SF_Tot_ep1 = new TH2D("h_SF_Tot_ep1", "", 200, 0., 0.95*DDVCSTools::Eb, 200, 0., 0.5);

    hists->h_SF_ECin_PCal_em1 = new TH2D("h_SF_ECin_PCal_em1", "", 200, 0., 0.35, 200, 0., 0.35);
    hists->h_SF_ECin_PCal_ep1 = new TH2D("h_SF_ECin_PCal_ep1", "", 200, 0., 0.35, 200, 0., 0.35);

    //------------------------------------ FS Histograms -----------------------------------
    hists->h_MXRecoil1 = new TH1D("h_MXRecoil1", "", 200, 0., 3.5);
    hists->h_MXRecoil2 = new TH1D("h_MXRecoil2", "", 200, 0., 3.5);
    hists->h_MXRecoil3 = new TH1D("h_MXRecoil3", "", 200, 0., 3.5);

    hists->h_vt_em1_1 = new TH1D("h_vt_em1_1", "", 200, -15., 15.);
    hists->h_vt_em2_1 = new TH1D("h_vt_em2_1", "", 200, -15., 15.);
    hists->h_vt_ep1 = new TH1D("h_vt_ep1", "", 200, -15., 15.);
    hists->h_vt_prot1 = new TH1D("h_vt_prot1", "", 200, -15., 15.);

    hists->h_Minv12_1 = new TH2D("h_Minv12_1", "", 200, 0., 3., 200, 0.,3.);
    hists->h_MxRecoil_Minv1_1 = new TH2D("h_MxRecoil_Minv1_1", "", 200, 0., 3., 200, 0.,2.);
    hists->h_MxRecoil_Minv2_1 = new TH2D("h_MxRecoil_Minv2_1", "", 200, 0., 3., 200, 0.,2.);
}

void InitCuts(AnaDB &db, int run, AnaCuts &cuts) {
    auto eTable = db.GetRunKeyedTable("DDVCS/ElectronDDVCS/AnaCuts");

    cuts.emVzMin = eTable->Get<double>(run, "emVzMin").value();
    cuts.emVzMax = eTable->Get<double>(run, "emVzMax").value();
    cuts.emVtMin = eTable->Get<double>(run, "emVtMin").value();
    cuts.emVtMax = eTable->Get<double>(run, "emVtMax").value();

    cuts.emChi2PIDMin = eTable->Get<double>(run, "emChi2PIDMin").value();
    cuts.emChi2PIDMax = eTable->Get<double>(run, "emChi2PIDMax").value();
    cuts.emPCalSFMin = eTable->Get<double>(run, "emPCalSFMin").value();
    cuts.emSFSumMin = eTable->Get<double>(run, "emSFSumMin").value();
    cuts.epVzMin = eTable->Get<double>(run, "epVzMin").value();
    cuts.epVzMax = eTable->Get<double>(run, "epVzMax").value();
    cuts.epVtMin = eTable->Get<double>(run, "epVtMin").value();
    cuts.epVtMax = eTable->Get<double>(run, "epVtMax").value();
    cuts.epPCalSFMin = eTable->Get<double>(run, "epPCalSFMin").value();
    cuts.epSFSumMin = eTable->Get<double>(run, "epSFSumMin").value();
    cuts.epChi2PIDMin = eTable->Get<double>(run, "epChi2PIDMin").value();
    cuts.epChi2PIDMax = eTable->Get<double>(run, "epChi2PIDMax").value();
    cuts.protVzMin = eTable->Get<double>(run, "protVzMin").value();
    cuts.protVzMax = eTable->Get<double>(run, "protVzMax").value();
    cuts.protVtMin = eTable->Get<double>(run, "protVtMin").value();
    cuts.protVtMax = eTable->Get<double>(run, "protVtMax").value();
    cuts.protChi2PIDMin = eTable->Get<double>(run, "protChi2PIDMin").value();
    cuts.protChi2PIDMax = eTable->Get<double>(run, "protChi2PIDMax").value();

    cuts.PCalEMin = eTable->Get<double>(run, "PCalEMin").value();
    cuts.PCalUmin = eTable->Get<double>(run, "PCalUmin").value();
    cuts.PCalVmin = eTable->Get<double>(run, "PCalVmin").value();
    cuts.PCalWmin = eTable->Get<double>(run, "PCalWmin").value();

    cuts.Eb = eTable->Get<double>(run, "Eb").value();
    DDVCSTools::Eb = cuts.Eb;
    std::cout<<"Kuku 1"<<std::endl;
}

bool ep_PCalECin_SF_Cut(const AnaCuts& cuts, const RecParticle& part) {
    return (part.energyECin() < 0.001 && part.SFPCal() > cuts.epPCalSFMin) || (part.energyECin() > 0.001 && (part.SFPCal() + part.SFECin() > cuts.epSFSumMin));
}


#endif //DDVCS_STUDIES_ANAELDDVCSHELPER_H
