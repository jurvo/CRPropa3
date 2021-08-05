#ifndef CRPROPA_DIFFUSIONCOEFFICENT_H
#define CRPROPA_DIFFUSIONCOEFFICENT_H

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/JF12Field.h"

#include <string>
#include <cmath>

namespace crpropa{

/**
 @class DiffusionTensor
 @brief Abstract base class for diffusion tensors
*/
class DiffusionTensor: public Referenced {
    private:
        std::string description;
    public:
        virtual ~DiffusionTensor(){
        }

        // diagonal entries of the diffusion tensor in the frame of the local magnetic fieldline. x direction is parallel to the fieldline, y and z are perpendicular
        virtual Vector3d getDiffusionKoefficent(Candidate *cand) const {
            return Vector3d(0.);
        };

        /*virtual double getKappaParallel(Candidate *cand){
            return 0;
        };
        
        virtual double getKappaPerpendicular(Candidate *cand){
            return 0;
        };
        
        virtual double getKappaPerpendicular2(Candidate *cand){
            return getKappaPerpendicular(cand);
        };*/
        virtual std::string getDescription() const{
            return "diffusion tensor";
        };
};


} // namespace


#endif // CRPROPA_DIFFUSIONCOEFFICENT_H