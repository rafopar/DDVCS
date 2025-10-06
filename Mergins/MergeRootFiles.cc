/* 
 * File:   MergeRootFiles.cc
 * Author: rafopar
 *
 * Created on February 26, 2025, 5:27 PM
 */

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
/* 
 * This code will read two input files that are created by the CLAS12Lund2Root, and will
 * make a new LUND file where every events of the LUND will contain one event from each input files.
 */

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if (argc < 5) {
        cout << "The command should look like " << endl;
        cout << "./MergeRootFiles.exe incRun GrapeRun NoEv NoEvPerFile"<<endl;
	cout<<" Exiting ..."<<endl;
	exit(1);
    }

    int incRun = atoi(argv[1]);
    int GrapeRun = atoi(argv[2]);
    int NeV = atoi(argv[3]);
    int NeVPerFile = atoi(argv[4]);

    TFile *file1 = TFile::Open(Form("InpData/Incluseive_Run_%d.root", incRun));
    TFile *file2 = TFile::Open(Form("InpData/grape_Run_%d.root", GrapeRun));

    if (!file1 || !file2 || file1->IsZombie() || file2->IsZombie()) {
        std::cerr << "Error: Could not open one or both ROOT files!" << std::endl;
        return 1;
    }

    // Access TTrees
    TTree *tree1 = (TTree*) file1->Get("tr1");
    TTree *tree2 = (TTree*) file2->Get("tr1");

    if (!tree1 || !tree2) {
        std::cerr << "Error: One or both TTrees not found!" << std::endl;
        return 1;
    }

    int iFile = 0;
    ofstream lund_Out(Form("OutData/IncRun_%d_GrapeRun_%d/merged_IncRun_%d_GrapeRun_%d_file_%d.txt", incRun, GrapeRun, incRun, GrapeRun, iFile));

    // Define variables to hold the data
    const int nMaxPart = 50;
    int nPart1;

    int A_Targ, Z_Targ, beamType, InterNuclID, ProcessID;
    double pol_targ, pol_beam, Eb, EvWeight;

    int index_1[nMaxPart];
    double t_live_1[nMaxPart]; // lifetime [ns]
    int type_1[nMaxPart]; // 1=active
    int pid_1[nMaxPart];
    int parentInd_1[nMaxPart];
    int daughtInd_1[nMaxPart];
    double px_1[nMaxPart];
    double py_1[nMaxPart];
    double pz_1[nMaxPart];
    double E_1[nMaxPart];
    double m_1[nMaxPart];
    double vx_1[nMaxPart];
    double vy_1[nMaxPart];
    double vz_1[nMaxPart];

    int nPart2;

    int index_2[nMaxPart];
    double t_live_2[nMaxPart]; // lifetime [ns]
    int type_2[nMaxPart]; // 1=active
    int pid_2[nMaxPart];
    int parentInd_2[nMaxPart];
    int daughtInd_2[nMaxPart];
    double px_2[nMaxPart];
    double py_2[nMaxPart];
    double pz_2[nMaxPart];
    double E_2[nMaxPart];
    double m_2[nMaxPart];
    double vx_2[nMaxPart];
    double vy_2[nMaxPart];
    double vz_2[nMaxPart];

    //    int A_Targ, Z_Targ, beamType, InterNuclID, ProcessID;
    //    double pol_targ, pol_beam, Eb, EvWeight;

    tree1->SetBranchAddress("A_Targ", &A_Targ);
    tree1->SetBranchAddress("Z_Targ", &Z_Targ);
    tree1->SetBranchAddress("pol_targ", &pol_targ);
    tree1->SetBranchAddress("pol_beam", &pol_beam);
    tree1->SetBranchAddress("beamType", &beamType);
    tree1->SetBranchAddress("Eb", &Eb);
    tree1->SetBranchAddress("InterNuclID", &InterNuclID);
    tree1->SetBranchAddress("ProcessID", &ProcessID);
    tree1->SetBranchAddress("EvWeight", &EvWeight);
    tree1->SetBranchAddress("nPart", &nPart1);
    tree1->SetBranchAddress("index", &index_1);
    tree1->SetBranchAddress("t_live", &t_live_1);
    tree1->SetBranchAddress("type", &type_1);
    tree1->SetBranchAddress("pid", &pid_1);
    tree1->SetBranchAddress("parentInd", &parentInd_1);
    tree1->SetBranchAddress("daughtInd", &daughtInd_1);
    tree1->SetBranchAddress("px", &px_1);
    tree1->SetBranchAddress("py", &py_1);
    tree1->SetBranchAddress("pz", &pz_1);
    tree1->SetBranchAddress("E", &E_1);
    tree1->SetBranchAddress("m", &m_1);
    tree1->SetBranchAddress("vx", &vx_1);
    tree1->SetBranchAddress("vy", &vy_1);
    tree1->SetBranchAddress("vz", &vz_1);

    tree2->SetBranchAddress("nPart", &nPart2);
    tree2->SetBranchAddress("index", &index_2);
    tree2->SetBranchAddress("t_live", &t_live_2);
    tree2->SetBranchAddress("type", &type_2);
    tree2->SetBranchAddress("pid", &pid_2);
    tree2->SetBranchAddress("parentInd", &parentInd_2);
    tree2->SetBranchAddress("daughtInd", &daughtInd_2);
    tree2->SetBranchAddress("px", &px_2);
    tree2->SetBranchAddress("py", &py_2);
    tree2->SetBranchAddress("pz", &pz_2);
    tree2->SetBranchAddress("E", &E_2);
    tree2->SetBranchAddress("m", &m_2);
    tree2->SetBranchAddress("vx", &vx_2);
    tree2->SetBranchAddress("vy", &vy_2);
    tree2->SetBranchAddress("vz", &vz_2);

    

    // Get total number of entries (assuming same length for both)
    int nEntries = std::min(tree1->GetEntries(), tree2->GetEntries());

    // Loop from 1 to 1M or up to the max available entries
    for (Long64_t i_ev = 0; i_ev < TMath::Min(NeV, nEntries ); i_ev++) {
    //for (Long64_t i_ev = 0; i_ev < nev; i_ev++) {

        if ((i_ev + 1) % NeVPerFile == 0) {

            lund_Out.close();
            iFile++;
            lund_Out.open(Form("OutData/IncRun_%d_GrapeRun_%d/merged_IncRun_%d_GrapeRun_%d_file_%d.txt", incRun, GrapeRun, incRun, GrapeRun, iFile));
        }


        tree1->GetEntry(i_ev);
        tree2->GetEntry(i_ev);

        //std::cout << "Entry " << i_ev << ": px1 " << px_1[0] << ", py1 " << py_1[0] << std::endl;

        int nPart = nPart1 + nPart2;

        lund_Out << nPart << setw(5) << A_Targ << setw(5) << Z_Targ << setw(5) << pol_targ << setw(9) << pol_beam << setw(5) << beamType << setw(10) << Eb << setw(9) << InterNuclID << setw(5) << ProcessID << setw(5) << EvWeight << endl;

        for (int ind = 0; ind < nPart1; ind++) {
            lund_Out << ind + 1 << setw(5) << t_live_1[ind] << setw(5) << type_1[ind] << setw(9) << pid_1[ind] << setw(5) << parentInd_1[ind] << setw(5) << daughtInd_1[ind] << setw(15)
                    << px_1[ind] << setw(15) << py_1[ind] << setw(15) << pz_1[ind] << setw(15) << E_1[ind] << setw(15) << m_1[ind] << setw(5) << vx_1[ind] << setw(5) << setw(5) << vy_1[ind] << setw(15) << vz_1[ind] << endl;
        }

        for (int ind = 0; ind < nPart2; ind++) {
            lund_Out << nPart1 + ind + 1 << setw(5) << t_live_2[ind] << setw(5) << type_2[ind] << setw(9) << pid_2[ind] << setw(5) << parentInd_2[ind] << setw(5) << daughtInd_2[ind] << setw(15)
                    << px_2[ind] << setw(15) << py_2[ind] << setw(15) << pz_2[ind] << setw(15) << E_2[ind] << setw(15) << m_2[ind] << setw(5) << vx_2[ind] << setw(5) << setw(5) << vy_2[ind] << setw(15) << vz_2[ind] << endl;
        }


    }

    // Close files
    file1->Close();
    file2->Close();

    delete file1;
    delete file2;



    return 0;
}

