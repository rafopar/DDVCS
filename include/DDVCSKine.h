//
// Created by rafopar on 10/30/25.
//

#ifndef DDVCSKINE_H
#define DDVCSKINE_H
#include <TLorentzVector.h>

#endif //DDVCSKINE_H

class DDVCSKine {
protected:
    DDVCSKine();

    virtual ~DDVCSKine();

    double fEb;
    double fMtarg;

    bool ftargetSet;
    bool fEbSet;
    bool fParticlesSet;
    bool fKinematicsComputed;

    double fMX2_Reaction;
    double fMx_Recoil;

    // Those LorentzVectors below are the same for both muon and electron DDVCS
    TLorentzVector *fL_beam, *fL_Targ, *fL_Recoil;

    TLorentzVector fL_MisReaction;                  // missing particle of the reaction (should be 0 in ideal case)

    void CheckKinematicsCalculated(const char* funcName) const;

private:

public:
    virtual void SetEb(double); // Sets the beam energy [GeV]
    virtual void SetMtarg(double); // Sets the mass of the target [GeV]

    double GetMxRecoil() const;

    static constexpr double fr2d = 57.2957795131;
};
