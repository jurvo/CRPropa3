#ifndef CRPROPA_PionDecay_H
#define CRPROPA_PionDecay_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Common.h"
#include "crpropa/Random.h"

using namespace crpropa;

class PionDecay: public Module
{
private:
    double mPion = 139.57039; // units MeV/c2
    double mMuon = 105.6583755; // units MeV/c2
    double r = pow_integer<2>(mMuon / mPion);
    double lambda = 1 - r;
    double h_nu_1, h_nue_1; // parameter for distribution of secondaries (following Kelner+ 2006)
    double samplingParameter = 1;

public:
    PionDecay();

    void process(Candidate* cand) const;
    void setSampling(double s);
    std::vector<double> sampleFraction() const;    
    void printSample();
    double f_e(double x) const;
    double f_nue(double x) const;
};

#endif // CRPROPA_PionDecay_H