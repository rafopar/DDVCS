

void NotWorkingCode(){

    int nPart;
    const int nMaxPart = 10;
    int index[nMaxPart];
    int pid[nMaxPart];
    int type[nMaxPart];
    int parentInd[nMaxPart];
    int daughtInd[nMaxPart];
    double t_live[nMaxPart];
    double px[nMaxPart];
    double py[nMaxPart];
    double pz[nMaxPart];
    double E[nMaxPart];
    double m[nMaxPart];
    double vx[nMaxPart];
    double vy[nMaxPart];
    double vz[nMaxPart];

    TFile *file_In = new TFile("tmp.root");
    TTree *tr1 = (TTree*)file_In->Get("tr1");

    tr1->SetBranchAddress("nPart", &nPart);
    tr1->SetBranchAddress("index", index);
    tr1->SetBranchAddress("t_live", t_live);
    tr1->SetBranchAddress("type", type);
    tr1->SetBranchAddress("pid", pid);
    tr1->SetBranchAddress("parentInd", parentInd);
    tr1->SetBranchAddress("daughtInd", daughtInd);
    tr1->SetBranchAddress("px", &px[0]);
    tr1->SetBranchAddress("py", py);
    tr1->SetBranchAddress("pz", pz);
    tr1->SetBranchAddress("E", E);
    tr1->SetBranchAddress("m", m);
    tr1->SetBranchAddress("vx", vx);
    tr1->SetBranchAddress("vy", vy);
    tr1->SetBranchAddress("vz", vz);

    int nEv = tr1->GetEntries();


        for( int iev = 0; iev < nEv; iev++ ){

        tr1->GetEntry(iev);

        cout<<"  ********************** Event "<<iev<<" **************************"<<endl;
        for( int i = 0; i < nPart; i++ ){
            cout<<index[i]<<"   "<<pid[i]<<"   "<<px[i]<<"   "<<py[i]<<"   "<<pz[i]<<"   "<<E[i]<<"   "<<m[i]<<"   "<<vx[i]<<"   "<<vy[i]<<"   "<<vz[i]<<endl;
        }

    }


}
