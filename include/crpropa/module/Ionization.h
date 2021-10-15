#ifndef CRPROPA_IONIZATION_H
#define CRPROPA_IONIZATION_H

#include "crpropa/Module.h"
#include "crpropa/massDistribution/Density.h"

using namespace crpropa;
/**
 * 
 * 
 * 
 */

class Ionization: public Module {

private:
    ref_ptr<Density> density;
    std::string propertyLabel = "ionizations";
    double IH = 13.598 * eV;
    double IH2 = 15.603 * eV;
    std::vector<double> anDiss = {-53.23133, 96.57789, -67.57069, 23.32707, -4.004618, 0.272652}; // parameter for polynomial fit
    std::vector<double> anDoub = {-125.8689, 172.0709, -108.8777, 34.18291, -5.358045, 0.335476}; // parameter for polynomial fit

public:
    Ionization(ref_ptr<Density> density);

    void process(Candidate *cand);


private: 
    /* crossections for the ionization after Padovani et al (A&A 501, 2, pp. 619-631, 2009)
        @param E: kinetic Energy of the incoming particle in [J]
        @return crossection for the process in [m^2]
    */
    double sigma_p_ion(double E);
    double sigma_p_ec(double E);
    double sigma_p_diss(double E);
    double sigma_p_doub(double E);

    double sigma_e_ion(double E);
    double sigma_e_diss(double E);
    double sigma_e_doub(double E);

    double energyLossRate(double E);
    void countInteractionsOnStep(double E, double step, Candidate *cand);

};






#endif // CRPROPA_IONIZATION_H