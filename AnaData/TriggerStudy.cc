//
// Created by rafopar on 10/15/25.
//

/*
 * The main objective of this code os to study events with electron trigger,
 * and check the fraction of events that don't have any track in the given sector.
 * This will help us to understand whether we should expect significantly higher trigger
 * rate compared to what we have estimated for DDVCS@uCLAS12 by selecting pions from the bit-8 data.
 */

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <iostream>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>


using namespace std;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: TriggerStudy.exe RUN" << std::endl;
        std::cerr << "Exiting" << std::endl;
        exit(1);
    }

    int run = atoi(argv[1]);
    char inputFile[256];

    const int nSec = 6;
    const int nMax_samePID = 15;
    const int HTCC_TYPE = 15;
    const int DET_DC = 6;

    sprintf(inputFile, "inpHipoFiles/rec_clas_%d.hipo", run);

    TFile *file_out = new TFile(Form("HistData/TriggerStudy_%d.root", run), "Recreate");

    TH1D h_trgBits1("h_trgBits1", "", 33, -0.5, 32.5);
    TH2D h_nTrk_Sec1("h_nTrk_Sec1", "", 11, -0.5, 10.5, nSec + 1, -0.5, nSec + 0.5);
    TH2D h_nNegTrk_Sec1("h_nNegTrk_Sec1", "", 11, -0.5, 10.5, nSec + 1, -0.5, nSec + 0.5);


    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;

    int evCounter = 0;

    hipo::bank bRunConf(factory.getSchema("RUN::config"));
    hipo::bank bRecPart(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));
    hipo::bank bRecTrk(factory.getSchema("REC::Track"));


    int ind_em[nMax_samePID];
    int ind_ep[nMax_samePID];
    int ind_prot[nMax_samePID];
    int sec_em[nMax_samePID];
    int sec_ep[nMax_samePID];


    try {
        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            if (evCounter % 10000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

            if (evCounter > 250000) {
                break;
            }

            event.getStructure(bRecPart);
            event.getStructure(bRecCalo);
            event.getStructure(bRecCC);
            event.getStructure(bRunConf);
            event.getStructure(bRecTrk);

            long trgWord = bRunConf.getLong("trigger", 0);
            int evNum = bRunConf.getInt("event", 0);

            bool em_trgBits[nSec] = {0, 0, 0, 0, 0, 0};
            int n_trk[nSec] = {0, 0, 0, 0, 0, 0};
            int n_neg_trk[nSec] = {0, 0, 0, 0, 0, 0};

            bool em_sectTrg = false;

            for (int ibit = 0; ibit < 64; ibit++) {
                if ((trgWord & (1 << ibit)) != 0) {
                    h_trgBits1.Fill(ibit);

                    if (ibit >= 1 && ibit <= nSec) {
                        em_trgBits[ibit - 1] = true;
                        em_sectTrg = true;
                    }
                }
            }


            if (!em_sectTrg) { continue; }


            map<int, int> ind_HTCC;
            map<int, int> ind_PCal;
            map<int, int> ind_ECin;
            map<int, int> ind_ECout;
            map<int, int> ind_FTOF;

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


            int nTrk = bRecTrk.getRows();

            //cout<<"TrgWord is "<<trgWord<<endl;
            for (int itrk = 0; itrk < nTrk; itrk++) {
                int status = bRecTrk.getInt("status", itrk);
                int detector = bRecTrk.getInt("detector", itrk);

                int q = bRecTrk.getInt("q", itrk);

                if (detector != DET_DC ) {
                    continue;
                }

                int sector = bRecTrk.getInt("sector", itrk) - 1;

              //  cout<<sector<<" ";
                n_trk[sector]++;

                if (q < 0) {
                    n_neg_trk[sector]++;
                }
            }
            //cout<<endl;

            for (int isec = 0; isec < nSec; isec++) {
                if ( em_trgBits[isec] ) {
                    h_nTrk_Sec1.Fill(n_trk[isec], isec);
                    h_nNegTrk_Sec1.Fill(n_neg_trk[isec], isec);

                    // if ( n_trk[isec] == 0 ) {
                    //     cout<<"Ev: "<<evNum<<"   Sec: "<<isec<<endl;
                    // }
                }
            }

        }
    } catch (const char *msg) {
        cerr << msg << endl;
    }


    gDirectory->Write();
    file_out->Close();

    return 0;
}
