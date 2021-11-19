#include "crpropa/module/Bremsstrahlung.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"


#include <cmath>

namespace crpropa
{

Bremsstrahlung::Bremsstrahlung(ref_ptr<Density> density, double limit, double threshold) : 
    density(density), limit(limit), havePhotons(false), secondaryThreshold(threshold), useDoubleIntegration(true), maximalVariationFactor(0.1) {
    tabDelta = {0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10};
    tabPhi1 = {45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65};
    tabPhi2 = {44.46, 44.38, 44.24, 43.65, 42.49, 40.19, 34.93, 29.78, 24.34, 17.28, 12.41};
    sigmaT = 6.652458558e-29 * meter * meter;
    alpha = 1/137.037;
    fMin = 1e-4;
}

void Bremsstrahlung::process(Candidate *cand) const {
    if(fabs(cand -> current.getId()) != 11)
        return; // only for electrons

    Vector3d pos = cand -> current.getPosition();
    double step = cand -> getCurrentStep();
    double n0 = getDensityAtPosition(pos);
    double Ein = cand -> current.getEnergy();

    // check for maximal variation of the density
    Vector3d dir = cand -> current.getDirection();
    Vector3d pos2 = pos + step*dir; 
    double n1 = getDensityAtPosition(pos2);
    if(abs(n1 - n0) / n0 > maximalVariationFactor) {
        // ToDo:    limit step size and perform half step 
        throw std::runtime_error("to big changes in density. Please take smaler step sizes. \n");
    }

    double crossection = getCrossection(Ein);

    if (crossection <= 0){
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
    double crossection = 0.;
    double epsMax = maximumEps(Ein);

    if(useDoubleIntegration) {
        double eps0 = epsBreak(Ein);
        crossection += gaussInt([this, Ein](double eps) {return this -> sigmaH(Ein, eps); }, 0, eps0);
        crossection += gaussInt([this, Ein](double eps) {return this -> sigmaH(Ein, eps); }, eps0 + 1e-5 * Ein, epsMax);
    }
    else{
        crossection += gaussInt([this, Ein](double eps) {return this -> sigmaH(Ein, eps); }, 0, epsMax);
    }
    return crossection;
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

void Bremsstrahlung::setUseDoubleIntegration(bool use) {
    this -> useDoubleIntegration = use;
}

void Bremsstrahlung::setMaximalVariationFactor(double factor) {
    this -> maximalVariationFactor = factor;
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

bool Bremsstrahlung::getUseDoubleIntegration() const {
    return this -> useDoubleIntegration;
}

double Bremsstrahlung::getMaximumVariationFactor() const {
    return this -> maximalVariationFactor;
}

double Bremsstrahlung::samplePhotonEnergy(double Ein) const {
    double epsMin = 0;
    double epsMax = maximumEps(Ein); 
    Random &random = Random::instance();

    for (int i = 0; i < 1000000; i++)
    {
        double eps = epsMin + random.rand() * (epsMax - epsMin);
        double pEps = sigmaH(Ein, eps);
        if(pEps > random.rand() * sigmaT/eV){
            return eps;
        }
    }
    throw std::runtime_error("error in sampling: no photon found! \n");    
}

// parametrisation of crossection after Schlickeiser 2002
double Bremsstrahlung::sigmaH(double E, double Egamma) const {
    double Delta = Egamma * mc2 / 4 / alpha / E / (E - Egamma);
    double Phi1, Phi2;
    if (Delta < 2) {
        Phi1 = interpolate(Delta, tabDelta, tabPhi1);
        Phi2 = interpolate(Delta, tabDelta, tabPhi2);
    }
    else {
        Phi1 = 4 * (log(2 * E / mc2 * (E - Egamma) / Egamma) - 1. / 2.);
        Phi2 = Phi1;
    }
    double sigma = (1 + pow_integer<2>(1 - Egamma / E) ) * Phi1 - 2. / 3. * (1 - Egamma / E) * Phi2;
    sigma *= 3. / 8. / M_PI * alpha * sigmaT;
    return sigma / Egamma / eV;
}

double Bremsstrahlung::maximumEps(double E) const {
    return E / (mc2 / 2 / E * exp(0.5) + 1);
}

double Bremsstrahlung::epsBreak(double E) const {
    return E / (mc2 / 8. / alpha / E + 1);
}


} // namespace crpropa
