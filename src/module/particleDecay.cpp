# include "crpropa/module/ParticleDecay.h"
#include "crpropa/Random.h"

ParticleDecay::ParticleDecay()
{
    // decay in pion rest frame
    E_mu_0 = (mPion * mPion + mMuon * mMuon) / 2 / mPion * c_squared;
    E_nu_0 = (mPion * mPion - mMuon * mMuon) / 2 / mPion * c_squared;
    p_mu_0 = std::sqrt(E_mu_0 * E_mu_0 - mMuon * mMuon * c_squared * c_squared) / c_light;
    p_nu_0 = E_nu_0 / c_light; // neglecting neutrino mass
}

void ParticleDecay::process(Candidate* cand) const {
    int id = cand -> current.getId();
    if (fabs(id) != 211)
        return; // only for pions 
    
    double gamma = cand -> current.getEnergy() / mPion; 
    double beta = std::sqrt(1 - gamma * gamma);
    Random random = Random::instance();
    double opticalDepth = std::exp(- cand -> getNextStep() / c_light / gamma  / lifetimePion);
    if (random.rand() < opticalDepth)
        return; // no decay
    
    // boosting decay products from pion rest frame (precalculated) into lab frame
    double theta = random.rand() * 2 * M_PI;
    double cTheta = std::cos(theta);

    double E_mu = gamma * E_mu_0 - gamma * beta * p_mu_0 * cTheta;
    double E_nu = gamma * E_nu_0 - gamma * beta * p_nu_0 * cTheta;

    // creating secondary neutrino
    int sign = id / fabs(id); 
    cand -> addSecondary(sign * 14, E_nu);
    
    // deactivate primary
    cand -> setActive(false);

    // add muon
    cand -> addSecondary(- sign * 13, E_nu);
}