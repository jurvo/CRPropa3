#ifndef CRPROPA_BREMSSTRAHLUNG_H
#define CRPROPA_BREMSSTRAHLUNG_H

#include "crpropa/Module.h"
#include "crpropa/Vector3.h"
#include "crpropa/massDistribution/Density.h"
#include <vector>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class Bremsstrahlung
 @brief Losses due Bremsstrahlung with gas.

 The differential crossection is parametrized after Schlickeiser 2002. It uses the scattering functions Phi_1 and Phi_2 after Blumenthal and Gould (Rev. Mod. Phys. 42, 237 (1970)). 
 For the highest photon energies the approximation of unschielded charge is used (Δ > 2).
 The density distribution considers only hydrogen, where the corssection for molecular hydrogen (H2) is 2*sigma_HI. The error of this approximation is less than 3%.
 */

class Bremsstrahlung: public Module{
protected:
    ref_ptr<Density> density; // density distribution to interact with
    double limit;   // limit next step length to limit the interaction probability to a given value
    bool havePhotons; // create secondary photons
    double secondaryThreshold; // threshold energy of the secondary photons

    double alpha; // fine structure constant
    double sigmaT; // thomson crossection
    double mc2 = mass_electron * c_squared; // rest energy of a electron

public:
    Bremsstrahlung(
        ref_ptr<Density> density, 
        double limit= 0.01, 
        double secondaryThreshold = 10 * keV
    );

    void process(Candidate *candidate) const;

    // differential crossection. Parametrized after Schlickeiser 2002.
    double differentialCrossection(double E, double eps) const;
    double getCrossection(double E) const;

    // calculate the produced photon energy. This is also used for the energy loss.
    double samplePhotonEnergy(double Ein) const;

    // set the density distribution
    void setDensity(ref_ptr<Density> density);
    ref_ptr<Density> getDensity() const;

    // set the limit of interaction probability per step
    void setLimit(double limit);
    double getLimit() const;

    // if true secondaries are produced
    void setHavePhotons(bool photons);
    bool getHavePhotons() const;

    // lower threshold for the produced secondaries. Requires that setHavePhotons is set to true.
    void setSecondaryThreshold(double threshold);
    double getSecondaryThreshold() const;

    // calculate density as n = n_HI + 2 * n_H2 
    double getDensityAtPosition(Vector3d &pos) const;
};

} // namespace crpropa

#endif // CRPROPA_BREMSSTRAHLUNG_H