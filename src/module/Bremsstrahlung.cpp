#include "crpropa/module/Bremsstrahlung.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"


#include <cmath>

namespace crpropa
{

Bremsstrahlung::Bremsstrahlung(ref_ptr<Density> density, double limit) : 
    density(density), limit(limit), havePhotons(false), secondaryThreshold(1*eV) {
    tabDelta = {0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10};
    tabPhi1 = {45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65};
    alpha = 1/137.037;
}

void Bremsstrahlung::process(Candidate *cand) const{
    if(fabs(cand -> current.getId()) != 11)
        return; // only for electrons

    // calculate energy loss following Schlickeiser 2002 eq. 4.4.19
    double lf = cand -> current.getLorentzFactor();
    double Delta = 1. / 4. /alpha / lf;
    
    double phi1 = interpolate(Delta, tabDelta, tabPhi1);
    Vector3d pos = cand -> current.getPosition();
    double nucleonDensity = density -> getNucleonDensity(pos);

    double dGdT =  3.9 * alpha * c_light * sigmaT / 8 / M_PI * lf * phi1 * nucleonDensity;
    double dG = dGdT * (cand -> getNextTimeStep());
    // apply energy loss and limit next step
    cand -> current.setLorentzFactor(lf - dG);
    cand -> limitNextTimeStep(limit* lf/dGdT);

    // optionally add secondary photons
    if(not(havePhotons))
        return;
}

void Bremsstrahlung::setDensity(ref_ptr<Density> density){
    this -> density = density;
}

void Bremsstrahlung::setLimit(double limit){
    this -> limit = limit;
}

void Bremsstrahlung::setHavePhotons(bool photons){
    this -> havePhotons = photons;
}

void Bremsstrahlung::setSecondaryThreshold(double threshold){
    this -> secondaryThreshold = threshold;
}

bool Bremsstrahlung::getHavePhotons() const{
    return this -> havePhotons;
}

double Bremsstrahlung::getSecondaryThreshold() const{
    return this -> secondaryThreshold;
}

double Bremsstrahlung::getLimit() const{
    return this -> limit;
}

} // namespace crpropa
