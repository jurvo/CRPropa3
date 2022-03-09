#include "crpropa/diffusionTensor/DiffusionTensor.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"

#include "kiss/logger.h"

using namespace crpropa;

void DiffusionTensorForParticles::add(ref_ptr<DiffusionTensor> tensor, int id) {
    tensorList.push_back(tensor);
    idList.push_back(id);
}

Vector3d DiffusionTensorForParticles::getDiffusionKoefficent(Candidate *cand) const {
    int id = cand -> current.getId();
    for(int i = 0; i < idList.size(); i++) {
        if(id == idList[i])
            return tensorList[i] -> getDiffusionKoefficent(cand);
    }

    KISS_LOG_WARNING << "DiffusionTensorForParticles: id not found in idList. Use defaultTensor instead. If defaultTensor is not specified the first tensor of the list will be used.\n";
    return Vector3d(0.); // TO DO:  Update with default tensor
}

std::string DiffusionTensorForParticles::getDescription() const {
    std::stringstream s;
    s << "Use DiffusionTensorForParticles with " << idList.size() << " tensors as listed below: \n";
    for(int i = 0; i < idList.size(); i++) {
        s << "For id: " << idList[i] << " use diffusionTensor: \n" << tensorList[i] -> getDescription() <<"\n\n"; 
    }
    return s.str();
}

// -----------------------------------------------------------------------------------------------------------

DiffusionTensorPowerlaw::DiffusionTensorPowerlaw(double kappa0, double rigidity, double alpha):
    kappa0(kappa0), rigidity(rigidity), alpha(alpha){}

Vector3d DiffusionTensorPowerlaw::getDiffusionKoefficent(Candidate *cand) const {
    double rig = cand->current.getEnergy() / cand->current.getCharge();
    double kappa = kappa0 * pow(rig/rigidity, alpha);
    return Vector3d(kappa);
}

double DiffusionTensorPowerlaw::getKappa0() const {
    return kappa0;
}

double DiffusionTensorPowerlaw::getRigidity() const {
    return rigidity;
}

double DiffusionTensorPowerlaw::getAlpha() const {
    return alpha;
}

void DiffusionTensorPowerlaw::setKappa0(double kappa0) {
    this -> kappa0 = kappa0;
}

void DiffusionTensorPowerlaw::setRigidity(double rig) {
    this -> rigidity = rig;
}

void DiffusionTensorPowerlaw::setAlpha(double alpha) {
    this -> alpha = alpha;
}

std::string DiffusionTensorPowerlaw::getDescription() const {
    std::stringstream s;
    s << "DiffusionTensorPowerlaw: Isotropic diffusion with a powerlaw.\n";
    s << "The diffusioncoefficent is normed to " << kappa0 << "m^2 / sec at rigidity of ";
    s << rigidity / mega / volt << " MV \n";
    s << "using a powerlaw index of " << alpha << "\n";
    return s.str();
}

// -----------------------------------------------------------------------------------------------------------

DiffusionTensorBrokenPowerlaw::DiffusionTensorBrokenPowerlaw(double kappa0, double rigRef, double rigBreak, double alphaLow, double alphaHigh) :
    kappa0(kappa0), rigidityReference(rigRef), rigidityBreak(rigBreak), alphaLow(alphaLow), alphaHigh(alphaHigh) 
{
    if(rigRef > rigBreak) {
        throw std::runtime_error("reference rigidity must be smaler than the break. \n");
    }
    kappaAtBreak = kappa0 * pow(rigBreak/rigRef, alphaLow);
}

Vector3d DiffusionTensorBrokenPowerlaw::getDiffusionKoefficent(Candidate *cand) const {
    double rig = cand -> current.getEnergy() / cand -> current.getCharge();
    double kappa;
    if(rig < rigidityBreak) {
        kappa = kappa0 * pow(rig / rigidityReference, alphaLow);
    }
    else{
        kappa = kappaAtBreak * pow(rig / rigidityBreak, alphaHigh);
    }
    return Vector3d(kappa);
}

void DiffusionTensorBrokenPowerlaw::setKappa0(double kappa){
    this->kappa0 = kappa;
}

void DiffusionTensorBrokenPowerlaw::setRigidityReference(double rig) {
    this->rigidityReference = rig;
}

void DiffusionTensorBrokenPowerlaw::setRigidityBreak(double rig) {
    this->rigidityBreak = rig;
}

void DiffusionTensorBrokenPowerlaw::setAlphaLow(double alpha) {
    this->alphaLow = alpha;
}

void DiffusionTensorBrokenPowerlaw::setAlphaHigh(double alpha) {
    this -> alphaHigh = alpha;
}

double DiffusionTensorBrokenPowerlaw::getKappa0() const {
    return kappa0;
}

double DiffusionTensorBrokenPowerlaw::getRigitityRefference() const {
    return rigidityReference;
}

double DiffusionTensorBrokenPowerlaw::getRigidityBreak() const {
    return rigidityBreak;
}

double DiffusionTensorBrokenPowerlaw::getAlphaLow() const {
    return alphaLow;
}

double DiffusionTensorBrokenPowerlaw::getAlphaHigh() const {
    return alphaHigh;
}

std::string DiffusionTensorBrokenPowerlaw::getDescription() const {
    std::stringstream s;
    s << "Diffusion tensor as a broken powerlaw. \n";
    s << "Using norm value of " << kappa0 << "m^2 / sec at rigidity of " << rigidityReference / mega / volt << "MV \n";
    s << "Break at " << rigidityBreak / mega / volt << " MV.\n";
    s << "Slope before break: " << alphaLow << "\t and after break: " << alphaHigh << "\n";
    return s.str(); 
}

