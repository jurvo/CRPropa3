#ifndef CRPROPA_IONIZATION_H
#define CRPROPA_IONIZATION_H

#include "crpropa/Module.h"
#include "crpropa/massDistribution/Density.h"

namespace crpropa {

class Ionization: public Module {

private:
    ref_ptr<Density> density;
    std::string propertyLabel = "ionizations";
    double limit;
    bool count;

    double IH = 13.598 * eV;
    double IH2 = 15.603 * eV;
    std::vector<double> anDiss = {-53.23133, 96.57789, -67.57069, 23.32707, -4.004618, 0.272652}; // parameter for polynomial fit
    std::vector<double> anDoub = {-125.8689, 172.0709, -108.8777, 34.18291, -5.358045, 0.335476}; // parameter for polynomial fit
    double sigmaT = 2 / M_PI / 3 * pow_integer<2>(fine_structure * h_planck / mass_electron / c_light);

public:
    Ionization(ref_ptr<Density> density, double limit = 0.01, bool count = true);

    void process(Candidate *cand) const;

    /* crossections for the ionization after Padovani et al (A&A 501, 2, pp. 619-631, 2009)
        @param E: kinetic Energy of the incoming particle in [J]
        @return crossection for the process in [m^2]
    */
    
    double sigmaPion(double E) const;
    double sigmaPec(double E) const;
    double sigmaPdiss(double E) const;
    double sigmaPdoub(double E) const;

    double sigmaEion(double E) const;
    double sigmaEdiss(double E) const;
    double sigmaEdoub(double E) const;

    double energyLossRateElectron(double E, Vector3d &pos) const;
    double energyLossRateProton(double E, Vector3d &pos) const; 
    void countInteractionsOnStep(double E, double step, Candidate *cand) const;
};

} // namespace

#endif // CRPROPA_IONIZATION_H