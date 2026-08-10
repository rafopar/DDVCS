//
// Created by rafopar on 6/4/26.
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

using namespace std;

int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Please provide the run number." << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    int run = atoi(argv[1]);
    char inputFile[256];

    sprintf(inputFile, "SkimDDVCS/hipo/elDDVCS_%06d.hipo", run);

    // Reading the DB
    AnaDB db("/work/clas12/rafopar/DB/RafoAnaDB.db");

    hipo::reader reader;
    reader.open(inputFile);
    hipo::dictionary dict;
    reader.readDictionary(dict);

    RecEvent recEvent(dict);
    hipo::event event;




}