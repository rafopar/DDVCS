//
// Created by rafopar on 10/30/25.
//

#ifndef ELECTRONDDVCSKINE_H
#define ELECTRONDDVCSKINE_H

#include <DDVCSKine.h>

class TLorentzVector;

class ElectronDDVCSKine : public DDVCSKine{

    public:
    ElectronDDVCSKine();

    /*
     * There are multiple ways of defining the kinematics, Either providing all 4 FS particle LorentzVectors or
     * providing at least three of them, and the missing one will be calculated using the four momentum conservation.
     */
    bool SetKin4FSParticles( TLorentzVector* L_em1, TLorentzVector* L_em2, TLorentzVector* L_ep, TLorentzVector* L_recoil );
    bool SetKineem1em2ep(TLorentzVector* L_em1, TLorentzVector* L_em2, TLorentzVector* L_ep);
    bool SetKineem1em2Recoil(TLorentzVector* L_em1, TLorentzVector* L_em2, TLorentzVector* L_Recoil);
    bool SetKineemepRecoil(TLorentzVector* L_em, TLorentzVector* L_ep, TLorentzVector* L_Recoil);

    /*
     * Getter methods
     */
    double GetMx2_Reaction() const;
    double GetQ2_1() const;
    double GetQ2_2() const;
    double GetQp2_1() const;
    double GetQp2_2() const;
    double GetNue_1() const;
    double GetNue_2() const;
    double GetxB_1() const;
    double GetxB_2() const;
    double GettM() const;
    double GetMinv_1() const;
    double GetMinv_2() const;
    double GetW_1() const;
    double GetW_2() const;
    double GetXi_1() const;
    double GetXi_2() const;
    double GetXiPrime_1() const;
    double GetXiPrime_2() const;
    double GetXX_GPD_1() const;
    double GetXX_GPD_2() const;
    double GetPhi_LH_1() const;
    double GetPhi_LH_2() const;



    private:
    /*
 * In the ep->e-e-e+p reaction there are 2 e-. Because of this an ambiguity arises
 * which one is the scattered beam electron and which one is a product of the decay of the
 * timelike photon (or one from the BH process)
 *
 * Because of this for each kinematic variable that is calculates using electron's kinematic,
 * we will have two clas members: e.g. fXX_1 and f_XX2. where for fXX_1 (fXX_2) the em1 (em2) is the beam
 * electron and em2 (em1) is the electron from the timelike photon decay.
 * the
 */
    double fQ2_1 = 0, fQ2_2 = 0;
    double fQp2_1 = 0, fQp2_2 = 0;
    double fnue_1 = 0, fnue_2 = 0;
    double fxB_1 = 0, fxB_2 = 0;
    double ftM = 0;  // The t_M will be calculated using proton's kinematics, and hence no ambiguity.
    double fMinv_1 = 0, fMinv_2 = 0;
    double fW_1 = 0, fW_2 = 0;
    double fxi_1 = 0, fxi_2 = 0;
    double fxiPrime_1 = 0, fxiPrime_2 = 0;
    double fXX_GPD_1 = 0, fXX_GPD_2 = 0;
    double fPhi_LH_1 = 0, fPhi_LH_2 = 0;
    double fMX_recoil = 0; // This is the missing mass of the (e-e-e+) system, should peak at proton mass.

    TLorentzVector *fL_em1, *fL_em2, *fL_ep;
    void ComputeKinematics();

};


#endif //ELECTRONDDVCSKINE_H
