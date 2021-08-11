#ifndef CRPROPA_BOHMDIFFUSION_H
#define CRPROPA_BOHMDIFFUSION_H

#include "crpropa/diffusionTensor/DiffusioNTensor.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/Units.h"

#include <string>

namespace crpropa{

class BohmDiffusion: public DiffusionTensor {
    private:
        double kappa0;
        double normRig;

    public:
        BohmDiffusion(double kappa0 = 6.1e24, double normRig = 4e9*volt);

        Vector3d getDiffusionKoefficent(Candidate *cand) const;

        double getKappa0() const;
        double getNormRig() const;

        void setKappa0(double kappa);
        void setNormRig(double rig);

        std::string getDescription() const;
};

} // namespace

#endif // CRPROPA_BOHMDIFFUSION_H