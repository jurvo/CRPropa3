#include "crpropa/diffusionTensor/BohmDiffusion.h"

using namespace crpropa;

BohmDiffusion::BohmDiffusion(double kappa0, double normRig):
    kappa0(kappa0), normRig(normRig) {}

Vector3d BohmDiffusion::getDiffusionKoefficent(Candidate *cand) const{
    double rig = cand -> current.getRigidity();
    Vector3d tens(kappa0*rig/normRig);
    return tens;
}

double BohmDiffusion::getKappa0() const{
    return kappa0;
}

double BohmDiffusion::getNormRig() const{
    return normRig;
}

void BohmDiffusion::setKappa0(double k){
    kappa0 = k;
}

void BohmDiffusion::setNormRig(double r){
    normRig = r;
}

std::string BohmDiffusion::getDescription() const{
    std::stringstream ss;
    ss  << "Diffusion tensor for Bohm diffusion: \n"
        << "scaling the diffusioncoefficent with D = k0 * rig/normRig and using \n"
        << "k0: \t" << kappa0 << " m^2 / s \n"
        << "normRig:\t" << normRig << " volt\n";

    return ss.str();
}