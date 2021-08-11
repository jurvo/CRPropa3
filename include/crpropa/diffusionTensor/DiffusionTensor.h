#ifndef CRPROPA_DIFFUSIONCOEFFICENT_H
#define CRPROPA_DIFFUSIONCOEFFICENT_H

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"

#include <string>


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


} // namespace


#endif // CRPROPA_DIFFUSIONCOEFFICENT_H