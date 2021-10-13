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
    std::vector<double> tabDelta;
    std::vector<double> tabPhi1;
    std::vector<double> tabPhi2;
    double mc2 = mass_electron * c_squared;

public:
    Bremsstrahlung(ref_ptr<Density> density, double limit= 0.01, double secondaryThreshold = 10 * keV);

    void process(Candidate *candidate) const;

    void setDensity(ref_ptr<Density> density);
    void setLimit(double limit);
    void setHavePhotons(bool photons);
    void setSecondaryThreshold(double threshold);
    void setUseDoubleIntegration(bool use);
    void setMaximalVariationFactor(double factor);

    ref_ptr<Density> getDensity() const;
    double getLimit() const;
    bool getHavePhotons() const;
    double getSecondaryThreshold() const;
    bool getUseDoubleIntegration() const;
    double getMaximumVariationFactor() const;

private:
    double samplePhotonEnergy(double Ein) const;
    double sigmaH(double E, double Egamma) const;
    double maximumEps(double E) const;
    double epsBreak(double E) const;
};

} // namespace crpropa

#endif // CRPROPA_BREMSSTRAHLUNG_H