#ifndef CRPROPA_QUASILINEARDIFFUSION_H
#define CRPROPA_QUASILINEARDIFFUSION_H

#include "crpropa/diffusionTensor/DiffusionTensor.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/RealisticJF12.h"
#include "crpropa/Units.h"

namespace crpropa{

/**
 @class QLTDiffusion
 @brief Diffusion Tensor for the quasi-linear theory.  
 */
class QLTDiffusion: public DiffusionTensor {
    private:
        double epsilon; // ratio between perpendicular and parallel diffusion coefficent
        double kappa0;  // norm value for the diffusioncoefficent at given rigidity
        double alpha;   // spectral index of the diffusion coefficent
        double normRig; // rigidity to norm the diffusionCoefficent

    public:
        QLTDiffusion(double epsilon = 0.1 , double kappa0 = 6.1e24, double alpha = (1./3.), double normRig = 4e9*volt);
        
        Vector3d getDiffusionKoefficent(Candidate *cand) const;

        void setEpsilon(double epsilon);
        void setKappa0(double kappa0);
        void setAlpha(double alpha);
        void setNormRig(double rig);
        //void setDescription();

        double getEpsilon() const;
        double getAlpha() const;
        double getKappa0() const;
        double getNormRig() const;

	    std::string getDescription() const;
};

class QLTTurbulent: public DiffusionTensor{
    private:
	    ref_ptr<MagneticField> backgroundField;
        ref_ptr<TurbulentField> turbulentField;
        ref_ptr<RealisticJF12Field> fullField;
        double kappa0;      // value to norm the diffusioncoefficent at a given rigidity
        double alphaPara;   // spectral index for the parallel component
        double alphaPerp;   // spectral index for the perpendicular component
        double normTurbulence;  // value to norm the turbulence (probably at earth)
        double normRig;     // rigidity to norm the diffusioncoefficent
        bool useFullModel;

    public:
        QLTTurbulent(ref_ptr<MagneticField> background, ref_ptr<TurbulentField> turbulent, double kappa0 = 6.1e24, double alphaPara=(1./3.), double alphaPerp=(1./3.), double normRig=4.0e9);
        QLTTurbulent(ref_ptr<RealisticJF12Field> fullField, double kappa0 = 6.1e24, double alphaPara = (1./3.), double alphaPerp = (1./3.), double normRig = 4e9);

        Vector3d getDiffusionKoefficent(Candidate *cand) const;

        double getKappa0() const;
        double getAlphaPara() const;
        double getAlphaPerp() const;
        double getNormTurbulence() const;
        double getNormRigidity() const;

        void setKappa0(double kappa0);
        void setAlphaPara(double alpha);
        void setAlphaPerp(double alpha);
        void setAlpha(double alpha);
        void setNormTurbulence(double eta);
        void setNormRigidity(double rig);
        void normToEarthPosition(Vector3d posEarth = Vector3d(-8.5*kpc, 0., 0.));

        std::string getDescription() const;
};

class QLTRigidity: public DiffusionTensor{
    private:
        ref_ptr<MagneticField> backgroundField;
        ref_ptr<TurbulentField> turbulentField;

        ref_ptr<RealisticJF12Field> field;
        double kappa0;
        double normEta;
        double normRho; // reduced rigidity to norm
        double alphaPara;
        double alphaPerp;
        double correlationLength;
        Vector3d normPos; // position where the diffusion coefficent is normed. default at earth Vector3d(-8.5*kpc, 0, 0)
        //double calculateLamorRadius(ParticleState &state) const;
        bool useFullModel;

    public:
        QLTRigidity(ref_ptr<MagneticField> magField, ref_ptr<TurbulentField> turbField, double kappa0=6.1e24, double alphaPara=(1./3.), double alphaPerp=(1./3.));
        QLTRigidity(ref_ptr<RealisticJF12Field> field, double kappa0=6.1e24, double alphaPara = (1./3.), double alphaPerp= (1./3.));
        
        Vector3d getDiffusionKoefficent(Candidate *cand) const;

        void setMagneticField(ref_ptr<MagneticField> field);
        void setTurbulentField(ref_ptr<TurbulentField> field);
        void setKappa0(double kappa0);
        void setAlphaPara(double aPara);
        void setAlphaPerp(double aPerp);
        void setAlpha(double alpha);
        void normToPosition(const Vector3d &pos);
        void setNormEta(double eta);
        void setNormRho(double rho);

        ref_ptr<MagneticField> getMagneticField();
        ref_ptr<TurbulentField> getTurbulentField();
        double getKappa0() const;
        double getAlphaPara() const;
        double getAlphaPerp() const;
        double getNormEta() const;
        double getNormRho() const;
        Vector3d getNormPos() const;
        std::string getDescription() const;
};

} // namespace

#endif // CRPROPA_QUASILINEARDIFFUSION_H