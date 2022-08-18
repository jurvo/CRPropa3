#include "crpropa/module/PionDecay.h"
#include "kiss/logger.h"


PionDecay::PionDecay() {
    h_nu_1 = f_e(r);
    h_nue_1 = f_nue(r);
}

void PionDecay::setSampling(double s) {
    samplingParameter = s;
}

std::vector<double> PionDecay::sampleFraction() const {
    Random &random = Random::instance();
    

    std::vector<double> x;
    for(int i = 0; i < 4; i++){
        x.push_back(0.);
    }
    // first muon neutrino
    x[0] = random.rand() * lambda;

    double f_e_max = 2.5 * samplingParameter;
    double f_nu_max = 2.5 * samplingParameter;

    double test1, test2, test3;
    // sample in single step x1, x2 and x3
    do{
        x[1] = random.rand();
        x[2] = random.rand();
        x[3] = 1 - x[0] - x[1] - x[2];
        if(x[3] <= 0){
            continue; // not allowed value
        }

        test1 = random.rand() * f_e_max;
        test2 = random.rand() * f_nu_max;
        test3 = random.rand() * f_e_max;
    } while ((f_e(x[1]) > test1) & (f_nue(x[2]) > test2) & (f_e(x[3]) > test3));
    return x;
}


void PionDecay::process(Candidate* cand) const {
    int id = cand -> current.getId();
    if (fabs(id) != 211){ 
        return; // only for pions
    }
    double E_pi = cand -> current.getEnergy();
    
    // std::vector<double> x = sampleFraction();
    // sampling energy fraction of secondaries

    Random &random = Random::instance();
    
    std::vector<double> x = sampleFraction();

    // creating secondaries
    int sign = id / fabs(id); 
    cand -> addSecondary(sign * 14, x[0] * E_pi);   // first muon neutrino
    cand -> addSecondary(- sign * 11, x[1] * E_pi); // elektron, positron
    cand -> addSecondary(sign * 12, x[2] * E_pi);   // elektron neutrino
    cand -> addSecondary(- sign * 14, x[3] * E_pi); // second muon neutrino 

    // deactivate primary
    cand -> setActive(false); 
    
}

double PionDecay::f_e(double x) const {
    double f;
    if (x >= r) {
        f = (3 - 2 * r) / (9 * pow_integer<2>(1 - r)) * (9 * x * x - 6 * std::log(x) - 4 * pow_integer<3>(x) - 5);
    }
    if (r > x) {
        f = h_nu_1 + (1 + 2 * r) * (r - x) / (9 * r * r) * (9 * (r + x) - 4 * (r * r + r * x + x * x));
    }
    return f;
}

double PionDecay::f_nue(double x) const {
    double f;
    double einsMinusR = 1 - r;
    double eMinusX = 1 - x;
    if (x >= r) { // see Erratum PRD 79, 039901(E) (2009)
        double g1 = 2. / 3. / einsMinusR / einsMinusR;
        double g2 = eMinusX * (6 * eMinusX * eMinusX + r * (5 + 5 * x - 4 * x * x));
        double g4 = 6 * r * std::log(x);
        return g1 * (g2 + g4);
    }
    if (r > x) {
        double hb1 = 2 * (r - x) / 3. / r / r;
        double hb2 = 7*r*r - 4*r*r*r + 7*x*r - 4*x*r*r - 2*x*x - 4*x*x*r;
        return h_nue_1 + hb1 * hb2;
    }
}

void PionDecay::printSample() {
    std::vector<double> x = sampleFraction();
    for(int i = 0; i< 4; i++){
        std::cout << x[i] << "\n";
    }
}