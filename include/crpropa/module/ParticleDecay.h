#ifndef CRPROPA_PARTICLEDECAY_H
#define CRPROPA_PARTICLEDECAY_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"


using namespace crpropa;

class ParticleDecay: public Module
{
private:
    double lifetimePion = 2.6e-8;
    double mPion = 139.57039 * MeV / c_squared;
    double mMuon = 105.6583755 * MeV / c_squared;
    double E_mu_0, E_nu_0; // energy after decay in pion rest frame
    double p_mu_0, p_nu_0; // momentum after decay in pion rest frame
public:
    ParticleDecay();

    void process(Candidate* cand) const;
};

#endif // CRPROPA_PARTICLEDECAY_H