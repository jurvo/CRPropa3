#include "crpropa/diffusionTensor/TurbulentDiffusion.h"
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

TurbulentDiffusion::TurbulentDiffusion(ref_ptr<MagneticField> backgroundField, ref_ptr<TurbulentField> turbulentField, bool useNormValue):
    backgroundField(backgroundField), turbulentField(turbulentField), useNormValue(useNormValue){
        loadData(getDataPath("TurbulenceData.txt"));
        useFullModel = false;
        correlationLength = turbulentField->getCorrelationLength();
        if(useNormValue){
            Vector3d normPos(-8.5*kpc, 0., 0.); // default earth position
            normToPosition(normPos);
        }    
}

TurbulentDiffusion::TurbulentDiffusion(ref_ptr<JF12Field> JF12, bool useNormValue):
    JF12(JF12), useNormValue(useNormValue){
        useFullModel = true;
        loadData(getDataPath("TurbulenceData.txt"));
        correlationLength = 108.80007902210151 * parsec;
        if(useNormValue){
            Vector3d normPos(-8.5*kpc, 0., 0.);
            normToPosition(normPos);
        }
}

void TurbulentDiffusion::loadData(std::string filename){
    // first clear all old Data
    turbulence.clear();
    alphaPara.clear();
    alphaPerp.clear();
    kappaPara.clear();
    kappaPerp.clear();

	std::ifstream infile(filename.c_str());

    if(!infile.good())
        throw std::runtime_error("TurbulentDiffusion diffusiontensor: could not open "+filename);
    
    int i = 0;
    while (infile.good()) {
        if (infile.peek() != '#') {
            double eta, aPar, aPer, kPar, kPer, err; // err is for not used error values
        
            infile >> eta >> aPar >> err;
            infile >> kPar >> err;
            infile >> aPer >> err;
            infile >> kPer >> err;

            turbulence.push_back(eta/std::sqrt(1+eta*eta)); // neue Turbulenzvariante   
            alphaPara.push_back(aPar);
            alphaPerp.push_back(aPer);
            kappaPara.push_back(kPar);
            kappaPerp.push_back(kPer);
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
        i++;
    }
    infile.close();
}

void TurbulentDiffusion::normToPosition(Vector3d &pos){
    useNormValue = true;
    normPosition = pos;
    kappa0 = 6.1e24;

    double b, B;
    if(useFullModel){
        b = JF12 -> getTurbulentStrength(pos);
        B = JF12 -> getRegularField(pos).getR();
    }
    else{
        b = turbulentField -> getField(pos).getR();
        B = backgroundField -> getField(pos).getR();
    }
    normTurbulence = b/std::sqrt(B*B+b*b);
    normRho = 4e9*volt/B/c_light/correlationLength;
}

double TurbulentDiffusion::getTurbulence(Vector3d &pos){
    double b, B;
    if(useFullModel){
        b = JF12 -> getTurbulentStrength(pos);
        B = JF12 -> getRegularField(pos).getR();
    }
    else{
        b = turbulentField -> getField(pos).getR();
        B = backgroundField -> getField(pos).getR();
    }
    return b/std::sqrt(B*B+b*b);
}

double TurbulentDiffusion::getAlphaPara(Vector3d &pos){
    double a = interpolate(getTurbulence(pos), turbulence, alphaPara);
    return a;
}

double TurbulentDiffusion::getAlphaPerp(Vector3d &pos){
    return interpolate(getTurbulence(pos), turbulence, alphaPerp);
}

double TurbulentDiffusion::getKappaParallel(Candidate *cand){
    Vector3d pos = cand -> current.getPosition();
    double rho;
    if(useFullModel){
        rho = calculateLamorRadius(cand->current, JF12)/correlationLength;
    }
    else{
        rho = calculateLamorRadius(cand->current, backgroundField)/correlationLength;
    }
    double eta = getTurbulence(pos);
    double kPar= 0;
    if(useNormValue){
        kPar = kappa0;
        kPar *= pow_integer<2>(normTurbulence/eta);
        kPar *= pow(rho/normRho, getAlphaPara(pos));
    }
    else{
        // full dependence from P.Reichherzer et al
        double k0 = interpolate(eta, turbulence, kappaPara);
        kPar = pow(10, k0 + getAlphaPara(pos)*log(rho));
    }
    return kPar;
}

double TurbulentDiffusion::getKappaPerpendicular(Candidate *cand){
    Vector3d pos = cand -> current.getPosition();
    double rho;
    if(useFullModel){
        rho = calculateLamorRadius(cand->current, JF12)/correlationLength;
    }
    else{
        rho = calculateLamorRadius(cand->current, backgroundField)/correlationLength;
    }
    double eta = getTurbulence(pos);
    double kPer;
    if(useNormValue){
        kPer = kappa0*pow_integer<2>(eta*normTurbulence)*pow(rho/normRho, getAlphaPerp(pos));
    }
    else{
        // full dependence from P.Reichherzer et al
        double k0 = interpolate(eta, turbulence, kappaPerp);
        kPer = pow(10, k0 + getAlphaPerp(pos)*log(rho));
    }
    return kPer;
}

double TurbulentDiffusion::getKappaPerpendicular2(Candidate *cand){
    return getKappaPerpendicular(cand);
}

bool TurbulentDiffusion::getUseNormValue() const{
    return useNormValue;
}

double TurbulentDiffusion::getNormTurbulence() const{
    if(!useNormValue){
        KISS_LOG_WARNING << "try to get normTurbulence in TurbulentDiffusion diffusion tensor while useNormValue is switched off \n" 
        << "The returned value will not be used. If you want to use it please turn it on with setUseNormValue(bool use)\n";
    }
    return normTurbulence;
}

double TurbulentDiffusion::getNormRho() const{
    if(!useNormValue){
        KISS_LOG_WARNING << "try to get normRho in TurbulentDiffusion diffusion tensor while useNormValue is switched off \n" 
        << "The returned value will not be used. If you want to use it please turn it on with setUseNormValue(bool use)\n";
    }
    return normRho;
}

double TurbulentDiffusion::getKappa0() const{
    if(!useNormValue){
        KISS_LOG_WARNING << "try to get kappa0 in TurbulentDiffusion diffusion tensor while useNormValue is switched off \n" 
        << "The returned value will not be used. If you want to use it please turn it on with setUseNormValue(bool use)\n";
    }
    return kappa0;
}

bool TurbulentDiffusion::getUseFullModel() const{
    return useFullModel;
}

double TurbulentDiffusion::getCorrelationLength() const{
    return correlationLength;
}

Vector3d TurbulentDiffusion::getNormPosition() const{
    return normPosition;
}

void TurbulentDiffusion::setUseNormValue(bool use){
    useNormValue = use;
}

void TurbulentDiffusion::setNormTurbulence(double turb){
    normTurbulence = turb;
}

void TurbulentDiffusion::setNormRho(double rho){
    normRho = rho;
}

void TurbulentDiffusion::setKappa0(double kappa){
    kappa0 = kappa;
}

// debug function
void TurbulentDiffusion::printData() {
    for(int i = 0; i<turbulence.size(); i++){
        std::cout << "t: \t" << turbulence[i] << " \n"
        << "aPar: \t" << alphaPara[i] << " \n"
        << "aPer: \t" << alphaPerp[i] << " \n"
        << "kPar: \t" << kappaPara[i] << " \n"
        << "kPer: \t" << kappaPerp[i] << " \n \n";
    }
}