//
// Created by rafopar on 10/30/25.
//

#ifndef DDVCSKINE_H
#define DDVCSKINE_H
#include <Math/LorentzVector.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>

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

    ROOT::Math::PxPyPzEVector fL_beam, fL_Targ, fL_Recoil;
    ROOT::Math::PxPyPzEVector fL_MisReaction;

    void CheckKinematicsCalculated(const char* funcName) const;

private:

public:
    virtual void SetEb(double); // Sets the beam energy [GeV]
    virtual void SetMtarg(double); // Sets the mass of the target [GeV]

    double GetMxRecoil() const;

    static constexpr double fr2d = 57.2957795131;
};

#endif //DDVCSKINE_H
