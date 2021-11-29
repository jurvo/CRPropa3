#include "crpropa/module/Ionization.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"

#include "kiss/logger.h"

using namespace crpropa;

Ionization::Ionization(ref_ptr<Density> density, double limit, bool count) : 
    density(density), limit(limit), count(count) { }

void Ionization::process(Candidate *cand) const {
    int id = abs(cand -> current.getId());
    if ((id != nucleusId(1,1)) && (id != 11))
        return; // only for electrons and for protons
    bool electron = (id == 11);
    double step = cand -> getCurrentStep();
    double E = cand -> current.getEnergy();
    Vector3d pos = cand -> current.getPosition();

    double dE = 0;
    if (electron){
        dE = energyLossRateElectron(E, pos) * step;
    }
    else {
        dE = energyLossRateProton(E, pos) * step;
    }

    if(dE == 0){
        KISS_LOG_WARNING << "step without energy loss \n";
    }
    else if (dE < 1 * keV)
    {
        KISS_LOG_DEBUG << "very low energyloss with step " << step/pc << "pc \n";
    }
    cand -> current.setEnergy(E - dE);
    cand -> limitNextStep(limit * E / dE * step);
    if(count)
        countInteractionsOnStep(E, step, cand);
}


double Ionization::sigmaPion(double E) const {
    // fiting parameter
    double A = 0.71; 
    double B = 1.63;
    double C = 0.51;
    double D = 1.24;
    
    double x = mass_electron * E / mass_proton / IH;

    double factor = 4 * M_PI * pow_integer<2>(bohr_radius);
    double sigmaL = factor * C * pow(x, D);
    double sigmaH = factor * (A * log(1 + x) + B) / x;

    return 1 / (1 / sigmaL + 1 / sigmaH);
}

double Ionization::sigmaPec(double E) const {
    double A = 1.044; 
    double B = 2.88;
    double C = 0.016; 
    double D = 0.136;
    double F = 5.86;
    int N = 2; // number of electrons
    
    double x = E / IH;

    double sigma = 4 * M_PI * pow_integer<2>(bohr_radius) * A * N * pow_integer<2>(IH / IH2 * x);
    sigma /= C + pow(x, B) + D * pow(x, F);
    return sigma;
}

double Ionization::sigmaPdiss(double E) const {
    return sigmaEdiss(mass_electron * E / mass_proton);
}

double Ionization::sigmaPdoub(double E) const {
    return sigmaEdoub(mass_electron * E / mass_proton);
}

double Ionization::sigmaEion(double E) const {
    double t = E / IH2;
    int N = 2; // number of electrons

    double n = 2.4;
    double A1 = 0.74;
    double A2 = 0.87;
    double A3 = -0.60;

    double F = (1 - pow(t, 1 - n)) / (n - 1) * pow(2 / (1 + t), n/2.) * (1 - pow(t, 1 - n / 2)) / (n - 2);
    double G = (A1 * log(t) + A2 + A3 / t) / t;

    double sigma = 4 * M_PI * pow_integer<2>(bohr_radius) * N * pow_integer<2>(IH / IH2) * F * G;
    return sigma;
}

double Ionization::sigmaEdiss(double E) const {
    if(E < 18.1 * eV) 
        return 0;
    
    double logE = log(E / eV);

    double logSigma = 0;
    for (int i = 0; i < 6; i++)
    {
        logSigma += anDiss[i] * pow(logE, i);
    }

    double sigma = exp(logSigma) * 1e-18 * cm * cm;
    return sigma;
}

double Ionization::sigmaEdoub(double E) const {
    if(E < 51 * eV) 
        return 0;

    double logE = log(E / eV);

    double logSigma = 0;
    for (int i = 0; i < 6; i++)
    {
        logSigma += anDoub[i] * pow(logE, i);
    }

    double sigma = exp(logSigma) * 1e-18 * cm * cm;
    return sigma;
}

double Ionization::energyLossRateElectron(double E, Vector3d &pos) const {
    double mc2 = mass_electron * c_squared;
    double gamma = E / mc2;
    double nHI = density -> getHIDensity(pos);
    double nH2 = density -> getH2Density(pos); 

    // loss in neutral atomic matter 
    // calculatet after Schlickeiser 2002 (eq 4.5.2)
    double dGdT = nHI * (log(gamma) + 2. / 3. * log(mc2 / IH));
    dGdT += nH2 * (log(gamma) + 2. / 3. * log(mc2 / IH2));
    dGdT *= 9. / 4. * c_light * sigmaT;
    return dGdT * mc2;
}

double Ionization::energyLossRateProton(double E, Vector3d &pos) const {
    double nHI = density -> getHIDensity(pos);
    double nH2 = density -> getH2Density(pos);
    double nHII = density -> getHIIDensity(pos);

    double mc2 = mass_electron * c_squared; // rest energy electron
    double Mc2 = mass_proton * c_squared; // rest energy proton
    double gamma = E/Mc2;
    double beta = sqrt(1 - 1 / gamma / gamma);
    double beta2 = beta * beta;

    // interactions in fully ionized plasma. (only proton)
    // calculatet after Schlickeiser 2002 (eq. 5.3.34) 
    // in the limit of relativistic protons and cold plasmas (W_e = 1)
    double dEdT_ion = 30 * c_light * sigmaT * mc2 / beta * nHII;

    // interactions in neutral matter
    double T = E - Mc2;
    double Qmax =  (T + 2 * Mc2) / (1 + pow_integer<2>((mass_proton + mass_electron)*c_light) / (2 * mass_electron * T) );
    double BHI = log((2 * mc2 * beta2 * Qmax) / (pow_integer<2>(IH / gamma))) - 2 * beta2;
    double BH2 = log((2 * mc2 * beta2 * Qmax) / (pow_integer<2>(IH2 / gamma))) - 2 * beta2;

    double dEdT_ato = 3 * c_light * sigmaT * mc2 / 4 / beta;
    dEdT_ato *= nHI * BHI + nH2 * BH2;

    return dEdT_ion + dEdT_ato;
}

void Ionization::countInteractionsOnStep(double E, double step, Candidate *cand) const {
    double totalCrossection;
    int id = cand->current.getId();
    bool electron = (fabs(id) == 11);
    Vector3d pos = cand->current.getPosition();

    if (electron) {
        totalCrossection += sigmaEion(E); 
        totalCrossection += sigmaEdiss(E);
        totalCrossection += sigmaEdoub(E);
    }
    else {
        int A = massNumber(id);
        int Z = chargeNumber(id);
        double eps = E / A;
        totalCrossection += sigmaPion(eps); 
        totalCrossection += sigmaPec(eps);
        totalCrossection += sigmaPdiss(eps);
        totalCrossection += sigmaPdoub(eps);
        totalCrossection *= pow_integer<2>(Z); // scaling for higher nuclei
    }

    double dens = density -> getH2Density(pos) + density -> getHIDensity(pos);
    double tau = totalCrossection * step * dens;
    int nrItteration = 1;
    while(tau > 0.01) {
        step /= 2;
        nrItteration *= 2;
        tau = totalCrossection * step * dens;
    }
    
    Random &random = Random::instance();
    int nrCount = 0;
    for (int i = 0; i < nrItteration; i++) {
        if(random.rand() > tau){
            nrCount++;
        }
    }

    if (cand -> hasProperty(propertyLabel)) {
        int counted = cand -> getProperty(propertyLabel);
        nrCount += counted;
    }
    cand -> setProperty(propertyLabel, nrCount); 
}