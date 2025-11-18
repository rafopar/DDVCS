//
// Created by rafopar on 10/30/25.
//

#include "DDVCSKine.h"

#include <cstdlib>
#include <iostream>
#include <ostream>

DDVCSKine::DDVCSKine(): ftargetSet(false), fEbSet(false), fParticlesSet(false), fL_beam(nullptr), fL_Targ(nullptr),
                        fL_Recoil(nullptr), fL_MisReaction() {
    fEb = 0;
    fMtarg = 0;
    fMX2_Reaction = 0;
    fMx_Recoil = 0;
    fKinematicsComputed = false;
}

void DDVCSKine::SetEb(double aEb) {
    if (aEb < 0) {
        std::cout << "Beam energy can not be negative. Exiting..." << std::endl;
        exit(1);
    }
    fEb = aEb;
    fL_beam = new TLorentzVector(0., 0.,fEb, fEb);
    fEbSet = true;
}

void DDVCSKine::SetMtarg(double aMtarg) {
    if (aMtarg < 0) {
        std::cout << "Mtarg can not be negative. Exiting..." << std::endl;
        exit(1);
    }

    fMtarg = aMtarg;

    fL_Targ = new TLorentzVector(0., 0.,0, fMtarg);

    ftargetSet = true;
}

double DDVCSKine::GetMxRecoil() const {
    CheckKinematicsCalculated(__func__);
    return fMx_Recoil;
}

void DDVCSKine::CheckKinematicsCalculated(const char *funcName) const{
    if (!fKinematicsComputed) {
        throw std::runtime_error(
            std::string("Error: \"") + funcName + "()\" called before setting the kinematics");
    }
}

DDVCSKine::~DDVCSKine() {
     delete fL_beam;
     delete fL_Targ;
     delete fL_Recoil;
}