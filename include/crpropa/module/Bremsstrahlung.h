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
private:
    ref_ptr<Density> density; // density distribution to interact with
    double limit;   // limit next step length to limit the interaction probability to a given value
    bool havePhotons; // create secondary photons
    double secondaryThreshold; // threshold energy of the secondary photons
    bool useDoubleIntegration; // calculate total crossection with two integrals instead of one
    double maximalVariationFactor; // maximal variation of the density between this step and the next one

    double alpha; // fine structure constant
    double sigmaT; // thomson crossection
    std::vector<double> tabDelta; // tabulated scattering functions 
    std::vector<double> tabPhi1;
    std::vector<double> tabPhi2;
    double mc2 = mass_electron * c_squared;

public:
    Bremsstrahlung(ref_ptr<Density> density, double limit= 0.01, double secondaryThreshold = 10 * keV);

    void process(Candidate *candidate) const;

    // set the density distribution
    void setDensity(ref_ptr<Density> density);
    // set the limit of interaction probability per step
    void setLimit(double limit);
    // if true secondaries are produced
    void setHavePhotons(bool photons);
    // lower threshold for the produced secondaries. Requires that setHavePhotons is set to true.
    void setSecondaryThreshold(double threshold);
    // calculate the total crossection with to integrals. (Default is true)
    void setUseDoubleIntegration(bool use);
    // allows a maximal variation of gas density for step length
    void setMaximalVariationFactor(double factor);

    ref_ptr<Density> getDensity() const;
    double getLimit() const;
    bool getHavePhotons() const;
    double getSecondaryThreshold() const;
    bool getUseDoubleIntegration() const;
    double getMaximumVariationFactor() const;

private:
    // calculate the produced photon energy. This is also used for the energy loss.
    double samplePhotonEnergy(double Ein) const;
    // differential crossection. Parametrized after Schlickeiser 2002.
    double sigmaH(double E, double Egamma) const;
    // maximal energy of secondary photon.
    double maximumEps(double E) const;
    // photon energy of the break between tabulated scattering functions and unshielded approximation.
    double epsBreak(double E) const;
};

} // namespace crpropa

#endif // CRPROPA_BREMSSTRAHLUNG_H