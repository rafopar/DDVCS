//
// Created by rafopar on 10/30/25.
//
#include <ElectronDDVCSKine.h>

ElectronDDVCSKine::ElectronDDVCSKine() : fL_em1(nullptr), fL_em2(nullptr), fL_ep(nullptr) {


}

bool ElectronDDVCSKine::SetKin4FSParticles(TLorentzVector* aL_em1, TLorentzVector* aL_em2, TLorentzVector* aL_ep, TLorentzVector* aL_recoil) {

    /*
     * At this point target mass and the beam energy should already be set
     */

    if ( !fEbSet || !ftargetSet ) {
        std::cerr<<"At this point the beam energy and the target mass should already be set. Exiting..."<<std::endl;
        exit(1);
    }

    fL_em1 = aL_em1;
    fL_em2 = aL_em2;
    fL_ep = aL_ep;
    fL_Recoil = aL_recoil;

    fParticlesSet = true;

    // All particles are set, so the reaction kinematics can ba calculated.
    ComputeKinematics();

    return true;
}

/*
 * This method assumes the full kinematics of the reaction ep->e-e+p is known, and it will calculate
 * all kinematic variables
 */
void ElectronDDVCSKine::ComputeKinematics() {

    fL_MisReaction = *fL_beam + *fL_Targ - *fL_em1 - *fL_em2 - *fL_ep - *fL_Recoil;

    /*
     * Calculate kinematic variables that assume the scattered electron the em1
     */
    const TLorentzVector fL_nu1 = *fL_em1 - *fL_beam;  // spacelike photon when em1 is the scattered electron
    const TLorentzVector fL_q1 = *fL_em2 + *fL_ep;    // timelike photon calculated when em1 is the scattered electron
    const TVector3 v3_beam_Eprime1 = fL_beam->Vect().Cross(fL_em1->Vect());
    const TVector3 v3_q_qprime1 = fL_nu1.Vect().Cross(fL_q1.Vect());

    fPhi_LH_1 = fL_em1->Vect().Dot(v3_q_qprime1) > 0 ? v3_beam_Eprime1.Angle(v3_q_qprime1) * fr2d : 360. - v3_beam_Eprime1.Angle(v3_q_qprime1) * fr2d;

    fQ2_1 = -fL_nu1.M2();
    fQp2_1 = fL_q1.M2();
    fMinv_1 = fL_q1.M();
    fnue_1 = -fL_nu1.E();
    fxB_1 = fQ2_1/(2*fMtarg*fnue_1);
    fxiPrime_1 = fxB_1/(2-fxB_1);
    fxi_1 = fxiPrime_1*(fQ2_1 + fQp2_1 )/fQ2_1;
    fXX_GPD_1 = 2*fxiPrime_1 - fxi_1;

    /*
     * Calculate kinematic variables that assume the scattered electron the em2
     */
    const TLorentzVector fL_nu2 = *fL_em2 - *fL_beam;  // spacelike photon when em2 is the scattered electron
    const TLorentzVector fL_q2 = *fL_em1 + *fL_ep;    // timelike photon calculated when em2 is the scattered electron
    const TVector3 v3_beam_Eprime2 = fL_beam->Vect().Cross(fL_em2->Vect());
    const TVector3 v3_q_qprime2 = fL_nu2.Vect().Cross(fL_q2.Vect());

    fPhi_LH_2 = fL_em2->Vect().Dot(v3_q_qprime2) > 0 ? v3_beam_Eprime2.Angle(v3_q_qprime2) * fr2d : 360. - v3_beam_Eprime2.Angle(v3_q_qprime2) * fr2d;

    fQ2_2 = -fL_nu2.M2();
    fQp2_2 = fL_q2.M2();
    fMinv_2 = fL_q2.M();
    fnue_2 = -fL_nu2.E();
    fxB_2 = fQ2_2/(2*fMtarg*fnue_2);
    fxiPrime_2 = fxB_2/(2-fxB_2);
    fxi_2 = fxiPrime_2*(fQ2_2 + fQp2_2 )/fQ2_2;
    fXX_GPD_2 = 2*fxiPrime_2 - fxi_2;

    /*
     * Variables that don't depend on the order of electrons
     */
    fMX2_Reaction = fL_MisReaction.M2();
    ftM = (*fL_Recoil - *fL_Targ).M2();

    fKinematicsComputed = true;
}

double ElectronDDVCSKine::GetMx2_Reaction() const {

    CheckKinematicsCalculated(__func__);

    return fMX2_Reaction;
}

double ElectronDDVCSKine::GetQ2_1() const {
    CheckKinematicsCalculated(__func__);
    return fQ2_1;
}
double ElectronDDVCSKine::GetQ2_2() const {
    CheckKinematicsCalculated(__func__);
    return fQ2_2;
}

double ElectronDDVCSKine::GetQp2_1() const {
    CheckKinematicsCalculated(__func__);
    return fQp2_1;
}
double ElectronDDVCSKine::GetQp2_2() const {
    CheckKinematicsCalculated(__func__);
    return fQp2_2;
}
double ElectronDDVCSKine::GetNue_1() const {
    CheckKinematicsCalculated(__func__);
    return fnue_1;
}
double ElectronDDVCSKine::GetNue_2() const {
    CheckKinematicsCalculated(__func__);
    return fnue_2;
}
double ElectronDDVCSKine::GetxB_1() const {
    CheckKinematicsCalculated(__func__);
    return fxB_1;
}
double ElectronDDVCSKine::GetxB_2() const {
    CheckKinematicsCalculated(__func__);
    return fxB_2;
}
double ElectronDDVCSKine::GettM() const {
    CheckKinematicsCalculated(__func__);
    return ftM;
}
double ElectronDDVCSKine::GetMinv_1() const {
    CheckKinematicsCalculated(__func__);
    return fMinv_1;
}
double ElectronDDVCSKine::GetMinv_2() const {
    CheckKinematicsCalculated(__func__);
    return fMinv_2;
}
double ElectronDDVCSKine::GetW_1() const {
    CheckKinematicsCalculated(__func__);
    return fW_1;
}
double ElectronDDVCSKine::GetW_2() const {
    CheckKinematicsCalculated(__func__);
    return fW_2;
}
double ElectronDDVCSKine::GetXi_1() const {
    CheckKinematicsCalculated(__func__);
    return fxi_1;
}
double ElectronDDVCSKine::GetXi_2() const {
    CheckKinematicsCalculated(__func__);
    return fxi_2;
}
double ElectronDDVCSKine::GetXiPrime_1() const {
    CheckKinematicsCalculated(__func__);
    return fxiPrime_1;
}
double ElectronDDVCSKine::GetXiPrime_2() const {
    CheckKinematicsCalculated(__func__);
    return fxiPrime_2;
}
double ElectronDDVCSKine::GetXX_GPD_1() const {
    CheckKinematicsCalculated(__func__);
    return fXX_GPD_1;
}
double ElectronDDVCSKine::GetXX_GPD_2() const {
    CheckKinematicsCalculated(__func__);
    return fXX_GPD_2;
}
double ElectronDDVCSKine::GetPhi_LH_1() const {
    CheckKinematicsCalculated(__func__);
    return fPhi_LH_1;
}
double ElectronDDVCSKine::GetPhi_LH_2() const {
    CheckKinematicsCalculated(__func__);
    return fPhi_LH_2;
}