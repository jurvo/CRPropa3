#include "crpropa/module/Bremsstrahlung.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"

#include <cmath>

namespace crpropa
{

Bremsstrahlung::Bremsstrahlung(ref_ptr<Density> density, double limit, double threshold) : 
    density(density), limit(limit), havePhotons(false), secondaryThreshold(threshold) {
    sigmaT = 6.652458558e-29 * meter * meter;
    alpha = 1/137.037;
}

void Bremsstrahlung::process(Candidate *cand) const {
    if(fabs(cand -> current.getId()) != 11)
        return; // only for electrons

    Vector3d pos = cand -> current.getPosition();
    double step = cand -> getCurrentStep();
    double n0 = getDensityAtPosition(pos);
    double Ein = cand -> current.getEnergy();
    double crossection = getCrossection(Ein);

    if (crossection < 0){
        throw std::runtime_error("no crossection left");
    }

    // limit next step
    cand -> limitNextStep(limit / crossection / n0); 

    // decide on interaction
    Random &random = Random::instance();
    double prob = n0 * crossection * step;
    if(random.rand() > prob)
        return; // no interaction
    
    // interaction will take place, choose right photon and apply energy loss
    double eps = samplePhotonEnergy(Ein);
    cand -> current.setEnergy(Ein - eps);
    if(havePhotons && (eps > secondaryThreshold))
        cand -> addSecondary(22, eps);
}

double Bremsstrahlung::getCrossection(double Ein) const {
    return gaussInt([this, Ein](double eps) {return this -> differentialCrossection(Ein, eps); }, 0, Ein);
}

void Bremsstrahlung::setDensity(ref_ptr<Density> density) {
    this -> density = density;
}

void Bremsstrahlung::setLimit(double limit) {
    this -> limit = limit;
}

void Bremsstrahlung::setHavePhotons(bool photons) {
    this -> havePhotons = photons;
}

void Bremsstrahlung::setSecondaryThreshold(double threshold) {
    this -> secondaryThreshold = threshold;
}

ref_ptr<Density> Bremsstrahlung::getDensity() const{
    return this -> density;
}

double Bremsstrahlung::getLimit() const {
    return this -> limit;
}

bool Bremsstrahlung::getHavePhotons() const {
    return this -> havePhotons;
}

double Bremsstrahlung::getSecondaryThreshold() const {
    return this -> secondaryThreshold;
}

double Bremsstrahlung::samplePhotonEnergy(double Ein) const {
    double epsMin = 0;
    double epsMax = Ein; 
    Random &random = Random::instance();

    for (int i = 0; i < 1000000; i++)
    {
        double eps = epsMin + random.rand() * (epsMax - epsMin);
        double pEps = differentialCrossection(Ein, eps);
        if(pEps > random.rand() * sigmaT * 1e-3){
            return eps;
        }
    }
    throw std::runtime_error("error in sampling: no photon found! \n");    
}

// parametrisation of crossection after Schlickeiser 2002 in the strong shielding limit
double Bremsstrahlung::differentialCrossection(double E, double eps) const {
    return 3 / (8 * M_PI * eps) * alpha * sigmaT * 45 * (4./3. - 4/3 * eps / E + pow_integer<2>(eps / E));
}

double Bremsstrahlung::getDensityAtPosition(Vector3d &pos) const {
    double n = density -> getHIDensity(pos) + 2 * density -> getH2Density(pos);
    return n;
}

} // namespace crpropa 