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
    ref_ptr<Density> density;
    double alpha= 1/137.037; // fine structure constant
    double sigmaT = 1/6 / M_PI * pow_integer<4>(eplus /c_light) * pow_integer<2>(1/epsilon0 /mass_electron);
    std::vector<double> tabDelta;
    std::vector<double> tabPhi1;
    double limit;
    bool havePhotons; // create secondary photons
    double secondaryThreshold; // threshold energy of the secondary photons

public:
    Bremsstrahlung(ref_ptr<Density> density, double limit= 0.01);

    void process(Candidate *cand);

    void setDensity(ref_ptr<Density> density);
    void setLimit(double limit);
    void setHavePhotons(bool photons);
    void setSecondaryThreshold(double threshold);

    bool getHavePhotons() const;
    double getSecondaryThreshold() const;
    double getLimit() const;
};

} // namespace crpropa

#endif // CRPROPA_BREMSSTRAHLUNG_H