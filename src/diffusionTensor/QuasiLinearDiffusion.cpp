#include "crpropa/diffusionTensor/QuasiLinearDiffusion.h"
#include "crpropa/Units.h"
#include "kiss/logger.h"
#include "crpropa/Common.h"

#include <sstream>
#include <fstream>

using namespace crpropa;


inline double calculateLamorRadius(ParticleState &state, ref_ptr<MagneticField> backgroundField) {
    double fieldStrength;
    Vector3d pos = state.getPosition();
    fieldStrength = backgroundField->getField(pos).getR();
    return state.getMomentum().getR()/std::abs(state.getCharge())/fieldStrength;
}
inline double calculateLamorRadius(ParticleState &state, double fieldStrength) {
    return state.getMomentum().getR()/std::abs(state.getCharge())/fieldStrength;
}


// QLTDiffusion ------------------------------------------------------------------
QLTDiffusion::QLTDiffusion(double epsilon , double kappa0, double alpha, double normRig ): 
    epsilon(epsilon), kappa0(kappa0), alpha(alpha), normRig(normRig) {}

void QLTDiffusion::setEpsilon(double eps){
    epsilon = eps;
}
void QLTDiffusion::setKappa0(double kap){
    kappa0 = kap;
}
void QLTDiffusion::setAlpha(double alph){
    alpha = alph;
}
void QLTDiffusion::setNormRig(double rig){
    normRig = rig;
}

double QLTDiffusion::getAlpha() const{
    return alpha;
}
double QLTDiffusion::getKappa0() const{
    return kappa0;
}
double QLTDiffusion::getEpsilon() const{
    return epsilon;
}
double QLTDiffusion::getNormRig() const{
    return normRig;
}

Vector3d QLTDiffusion::getDiffusionKoefficent(Candidate *cand) const{
    double rig = cand -> current.getRigidity();

    Vector3d tens(1, epsilon, epsilon);
    tens *= kappa0 * pow(std::abs(rig)/normRig, alpha);
    return tens;
}

std::string QLTDiffusion::getDescription() const{
    std::stringstream s;
    s << "Diffusion tensor for quasi-linear-theory (QLT)  with: \n"
    << "epsilon: \t" << epsilon <<"\n"
    << "kappa0: \t"<< kappa0 << " m^2/s \n"
    << "alpha:  \t" << alpha << "\n"
    << "normRig: \t" << normRig << "\n";
    return s.str();
}

// QLTTurbulent ------------------------------------------------------------------

QLTTurbulent::QLTTurbulent(ref_ptr<MagneticField> backgroundField, ref_ptr<TurbulentField> turbulentField, double kappa0, double alphaPara, double alphaPerp, double normRig):
    backgroundField(backgroundField), turbulentField(turbulentField), kappa0(kappa0), alphaPara(alphaPara), alphaPerp(alphaPerp), normRig(normRig) {
    useFullModel = false;
    normToEarthPosition(); 
}

QLTTurbulent::QLTTurbulent(ref_ptr<RealisticJF12Field> fullField, double kappa0, double alphaPara, double alphaPerp, double normRig):
    fullField(fullField), kappa0(kappa0), alphaPara(alphaPara), alphaPerp(alphaPerp), normRig(normRig) {
    useFullModel = true;
    normToEarthPosition();
}

Vector3d QLTTurbulent::getDiffusionKoefficent(Candidate *cand) const    {
    double rig = cand -> current.getRigidity();
    Vector3d pos = cand -> current.getPosition();
    double eta;
    if(useFullModel){
        eta = fullField -> getTurbulenceOverRegular(pos);
    }
    else{
        double b = turbulentField -> getField(pos).getR();
        double B = backgroundField -> getField(pos).getR();
        // double eta = b/std::sqrt(b*b + B*B);  // new approach by Julien
        eta = b/B; // std approach (e.g. Partrick)
    }

    Vector3d tens(kappa0);
    tens.x *= pow_integer<2>(normTurbulence/eta)*pow(rig/normRig, alphaPara);
    tens.y *= pow_integer<2>(normTurbulence*eta)*pow(rig/normRig, alphaPerp);
    tens.z = tens.y;

    return tens;
}

double QLTTurbulent::getKappa0() const{
    return kappa0;
}
double QLTTurbulent::getAlphaPara() const{
    return alphaPara;
}
double QLTTurbulent::getAlphaPerp() const{
    return alphaPerp;
}
double QLTTurbulent::getNormTurbulence() const{
    return normTurbulence;
}
double QLTTurbulent::getNormRigidity() const{
    return normRig;
}

void QLTTurbulent::setKappa0(double kap){
    kappa0 = kap;
}
void QLTTurbulent::setAlphaPara(double a){
    alphaPara = a;
}
void QLTTurbulent::setAlphaPerp(double a){
    alphaPerp = a;
}
void QLTTurbulent::setAlpha(double a){
    alphaPara = a;
    alphaPerp = a;
}
void QLTTurbulent::setNormTurbulence(double eta){
    normTurbulence = eta;
}
void QLTTurbulent::setNormRigidity(double rig){
    normRig = rig;
}

void QLTTurbulent::normToEarthPosition(Vector3d posEarth){
    if(useFullModel){
        setNormTurbulence(fullField -> getTurbulenceOverRegular(posEarth));
    }
    else{
        double b = turbulentField -> getField(posEarth).getR();
        double B = backgroundField ->getField(posEarth).getR();
        setNormTurbulence(b/B); 
    }

}

std::string QLTTurbulent::getDescription() const{
    std::stringstream ss;
    ss  << "Diffusion Tensor for the quasi-linear-theory (QLT) with turbulence dependence \n \n"
        << "Energyscaling with spectral index:\n"
        << "alpha parallel \t" << alphaPara << "\n"
        << "alpha perpendicular \t" << alphaPerp << "\n \n"
        << "norming the turbulence at value: \t" << normTurbulence << "\n"
        << "and at rigidity: \t" << normRig <<" volt \n"
        << "norming the diffusion coefficent to: \t" << kappa0 << " m2/s \n"; 
    return ss.str();
}

// QLTRigidity ------------------------------------------------------------------

QLTRigidity::QLTRigidity(ref_ptr<MagneticField> magField, ref_ptr<TurbulentField> turbField, double kappa0, double alphaPara, double alphaPerp)
    : backgroundField(magField), kappa0(kappa0), alphaPara(alphaPara), alphaPerp(alphaPerp){
        setTurbulentField(turbField);
        useFullModel = false;
}

QLTRigidity::QLTRigidity(ref_ptr<RealisticJF12Field> field, double kappa0, double alphaPara, double alphaPerp)
    : field(field), kappa0(kappa0), alphaPara(alphaPara), alphaPerp(alphaPerp) {
    useFullModel = true;
    normToPosition(normPos);
}

void QLTRigidity::setMagneticField(ref_ptr<MagneticField> field){
    backgroundField = field;
    normToPosition(normPos);
    useFullModel = false;
}

void QLTRigidity::setTurbulentField(ref_ptr<TurbulentField> field){
    turbulentField = field;
    correlationLength = field -> getCorrelationLength();
    normToPosition(normPos);
    useFullModel = false;
}

void QLTRigidity::setKappa0(double kap){
    kappa0 = kap;
}
void QLTRigidity::setAlphaPara(double a){
    alphaPara = a;
}
void QLTRigidity::setAlphaPerp(double a){
    alphaPerp = a;
}
void QLTRigidity::setAlpha(double a){
    alphaPara = a;
    alphaPerp = a;
}

void QLTRigidity::normToPosition(const Vector3d &pos){
    normPos = pos;
    Vector3d B;
    if(useFullModel){
        normEta = field->getTurbulenceOverRegular(pos);
        correlationLength = 60 * parsec;
        B = field->getRegularField(pos);
    }
    else{
        Vector3d b = turbulentField -> getField(pos);
        B = backgroundField -> getField(pos);
        normEta = b.getR()/B.getR();
    }
    normRho = 4e9*volt/B.getR()/c_light/correlationLength;
}

void QLTRigidity::setNormEta(double eta){
    normEta = eta;
}

void QLTRigidity::setNormRho(double rho){
    normRho = rho;
}

ref_ptr<MagneticField> QLTRigidity::getMagneticField(){
    return backgroundField;
}

ref_ptr<TurbulentField> QLTRigidity::getTurbulentField(){
    return turbulentField;
}

double QLTRigidity::getKappa0() const{
    return kappa0;
}

double QLTRigidity::getAlphaPara() const{
    return alphaPara;
}

double QLTRigidity::getAlphaPerp() const{
    return alphaPerp;
}

double QLTRigidity::getNormEta() const{
    return normEta;
}

double QLTRigidity::getNormRho() const{
    return normRho;
}

Vector3d QLTRigidity::getNormPos() const{
    return normPos;
}

Vector3d QLTRigidity::getDiffusionKoefficent(Candidate *cand) const {
    Vector3d pos = cand -> current.getPosition();
    double rho, eta;
    if(useFullModel){
        rho = calculateLamorRadius(cand -> current, field->getRegularField(pos).getR());
        eta = field->getTurbulenceOverRegular(pos);
    }
    else{
        rho = calculateLamorRadius(cand -> current, backgroundField)/correlationLength;
        double b = turbulentField -> getField(pos).getR();
        double B = backgroundField -> getField(pos).getR();
        eta = b/B;
    }
    
    Vector3d tens(kappa0);
    tens.x *= pow_integer<2>(normEta/eta)*pow(rho/normRho, getAlphaPara());  // parallel component
    tens.y *= pow_integer<2>(normEta*eta)*pow(rho/normRho, getAlphaPerp());  // perpendicular component
    tens.z = tens.y; // perpendicular components are degenerated

    return tens;
}