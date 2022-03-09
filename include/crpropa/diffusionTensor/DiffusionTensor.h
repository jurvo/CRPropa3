#ifndef CRPROPA_DIFFUSIONCOEFFICENT_H
#define CRPROPA_DIFFUSIONCOEFFICENT_H

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"

#include <string>
#include <vector>


namespace crpropa{

/**
 @class DiffusionTensor
 @brief Abstract base class for diffusion tensors
*/
class DiffusionTensor: public Referenced {
    public:
        virtual ~DiffusionTensor(){
        }

        // diagonal entries of the diffusion tensor in the frame of the local magnetic fieldline. x direction is parallel to the fieldline, y and z are perpendicular
        virtual Vector3d getDiffusionKoefficent(Candidate *cand) const {
            return Vector3d(0.);
        };

        virtual std::string getDescription() const{
            return "diffusion tensor";
        };
};

class DiffusionTensorForParticles: public DiffusionTensor {
private:
    std::vector<int> idList;
    std::vector< ref_ptr<DiffusionTensor> > tensorList;
    ref_ptr<DiffusionTensor> defaultTensor = NULL;

public:
    void add(ref_ptr<DiffusionTensor> tensor, int id);
    Vector3d getDiffusionKoefficent(Candidate *cand) const;
    std::string getDescription() const;
};

class DiffusionTensorPowerlaw: public DiffusionTensor {
protected:
    double kappa0;
    double rigidity;
    double alpha;

public:
    DiffusionTensorPowerlaw(double kappa, double rigidity, double alpha);
    Vector3d getDiffusionKoefficent(Candidate *cand) const;

    double getKappa0() const;
    double getRigidity() const;
    double getAlpha() const;

    void setKappa0(double kappa);
    void setRigidity(double rigidity);
    void setAlpha(double alpha);

    std::string getDescription() const;
};


class DiffusionTensorBrokenPowerlaw: public DiffusionTensor {
private:
    double kappa0;
    double rigidityReference;
    double rigidityBreak;
    double alphaLow;
    double alphaHigh;

    double kappaAtBreak;
public:

    DiffusionTensorBrokenPowerlaw(double kappa0, double rigRef, double rigBreak, double alpha1, double alpha2);
    Vector3d getDiffusionKoefficent(Candidate *cand) const;
    std::string getDescription() const;

    void setKappa0(double kappa0);
    void setRigidityReference(double rigRef);
    void setRigidityBreak(double rigBreak);
    void setAlphaLow(double a);
    void setAlphaHigh(double a);

    double getKappa0() const;
    double getRigitityRefference() const;
    double getRigidityBreak() const;
    double getAlphaLow() const;
    double getAlphaHigh() const;
};

/*
class DiffusionTensorMultipleBreak: public DiffusionTensor {
private:
    double kappa0; // value to norm
    double alpha0; // first slope
    double rig0; // norm of kappa
    std::vector<double> kappa; // kappa value for break point
    std::vector<double> rig; // break positions
    std::vector<double> alpha; // slope after break

public:
    DiffusionTensorMultipleBreak(double kappa0, double alpha0, double rig0);
    void addBreak(double rig, double alpha);
};
*/

} // namespace


#endif // CRPROPA_DIFFUSIONCOEFFICENT_H