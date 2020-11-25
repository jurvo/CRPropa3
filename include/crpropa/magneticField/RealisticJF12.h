#ifndef CRPROPA_REALISTICJF12_H
#define CRPROPA_REALISTICJF12_H

#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/JF12FieldSolenoidal.h"
#include "crpropa/magneticField/CMZField.h"
#include "kiss/logger.h"

namespace crpropa {

class RealisticJF12Field: public MagneticField {
private:
    double lc = 54*pc;
    double reduction = 10.;
    

public:
    JF12FieldSolenoidal* JF12;
    CMZField* CMZ;

    RealisticJF12Field(){
        JF12 = new JF12FieldSolenoidal(3 * kpc, 0.5*kpc);
        JF12 -> deactivateOuterTransition();
        JF12 -> randomTurbulent(0);
        JF12 -> randomStriated(0);
        CMZ = new CMZField();
    };

    Vector3d getField(const Vector3d pos) const {
        Vector3d b(0.);
        b += JF12->getRegularField(pos);
        b += JF12->getTurbulentField(pos)/reduction;
        b += CMZ->getField(pos);
        return b;
    };

    Vector3d getRegularField(const Vector3d pos) const {
        Vector3d b(0.);
        b += JF12->getRegularField(pos);
        b += CMZ->getField(pos);
        return b;
    };

    Vector3d getTurbulentField(const Vector3d pos) const {
        return JF12 -> getTurbulentField(pos)/reduction;
    };

    void setReduction(double red) {
        reduction = red;
    };

    double getReduction() const {
        return reduction;
    };

}; // class RealisticJF12Field

} // namespace crpropa

#endif //CRPROPA_REALISTICJF12_H
