#ifndef CRPROPA_TURBULENTDIFFUSION_H
#define CRPROPA_TURBULENTDIFFUSION_H

#include "crpropa/diffusionTensor/DiffusionTensor.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/JF12Field.h"
#include "crpropa/magneticField/RealisticJF12.h"

#include <cmath>
#include <string>

namespace crpropa {

class TurbulentDiffusion: public DiffusionTensor{
    protected:
        // results by P.Reichherzer (MNRAS)
        std::vector<double> turbulence;
        std::vector<double> alphaPara;
        std::vector<double> alphaPerp;
        std::vector<double> kappaPara;  // data in log-space
        std::vector<double> kappaPerp;  // data in log-space

        // alternative norming
        bool useNormValue;
        double normTurbulence;
        double normRho;
        double kappa0;
        double correlationLength;
        Vector3d normPosition;

        // fields for local turbulence
        ref_ptr<MagneticField> backgroundField;
        ref_ptr<TurbulentField> turbulentField;
        bool useFullModel; // if true all turbulence values are calculated by JF12 field
        ref_ptr<RealisticJF12Field> JF12;

    public:
        TurbulentDiffusion(ref_ptr<MagneticField> background, ref_ptr<TurbulentField> turbulent, bool useNormValue=false);
        TurbulentDiffusion(ref_ptr<RealisticJF12Field> JF12, bool useNormValue=false);

        void loadData(std::string filename);
        void normToPosition(Vector3d &pos);

        Vector3d getDiffusionKoefficent(Candidate *cand) const;

        double getAlphaPara(Vector3d &pos) const;
        double getAlphaPerp(Vector3d &pos) const;
        double getTurbulence(Vector3d &pos) const;

        double getKappaParallel(Candidate *cand);
        double getKappaPerpendicular(Candidate *cand);
        double getKappaPerpendicular2(Candidate *cand);

        bool getUseNormValue() const;
        bool getUseFullModel() const;
        double getNormTurbulence() const;
        double getNormRho() const;
        double getKappa0() const;
        double getCorrelationLength() const;
        Vector3d getNormPosition() const;
        
        void setUseNormValue(bool use);
        void setNormTurbulence(double turb);
        void setNormRho(double rho);
        void setKappa0(double kappa);

        void printData();
};


class TurbulentCSR: public TurbulentDiffusion{
    public:
        Vector3d getDiffusionKoefficent(Candidate *cand) const;
};


}; // namespace
#endif // CRPROPA_TURBULENTDIFFUSION_H