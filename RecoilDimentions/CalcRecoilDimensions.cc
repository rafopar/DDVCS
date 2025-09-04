/* 
 * File:   CalcRecoilDimentions.cc
 * Author: rafopar
 *
 * Created on March 18, 2025, 11:26 PM
 */

#include <cstdlib>

using namespace std;

/*
 * 
 */
void CalcRecoilDimensions() {

    const int nLayer = 6;
    const int nSec = 3;
    
    const double minRadius = 15; // cm
    const double delta_R = 1; // cm
    const double th_Min = 40*TMath::DegToRad();
    const double th_Max = 70*TMath::DegToRad();
    const double pitch_C_Strips = 0.05; // cm
    const double pitch_Z_Strips = 0.05; // cm
    const double l_targ = 7.5; // cm
    
    int n_Tot_Strips = 0;
    double TotArea = 0;
    for( int il = 0; il < nLayer; il++ ){
        
        double cur_R = minRadius + il*delta_R;
        double length = l_targ + cur_R*( 1./tan(th_Min) - 1./tan(th_Max) );
        
        double perimeter = 2*TMath::Pi()*cur_R;
        
        double area = length*perimeter;
        TotArea = TotArea + area;
        int n_Z_Strips = int(perimeter/pitch_Z_Strips);
        
        int n_C_Strips = nSec*int(length/pitch_C_Strips);
        
        int n_Tot_StripsInLayer = n_C_Strips + n_Z_Strips;
        n_Tot_Strips = n_Tot_Strips + n_Tot_StripsInLayer;
        
        cout<<"******** Layer "<<il<<" ********"<<endl;
        cout<<"* "<<endl;
        cout<<"* Length "<<length<<" cm,  Perimeter "<<perimeter<<" cm"<<"    Area = "<<length*perimeter<<" cm2"<<endl;
        cout<<"* C Strips "<<n_C_Strips<<",  Z strips "<<n_Z_Strips<<"  All Strips "<<n_Tot_StripsInLayer<<endl;
        cout<<"* "<<endl;
        
    }

    cout<<"* "<<endl;
    cout<<"* Total number of strips is  "<<n_Tot_Strips<<endl;
    cout<<"* Total Area is  "<<TotArea<<" cm2"<<endl;
    cout<<"* "<<endl;
    
}