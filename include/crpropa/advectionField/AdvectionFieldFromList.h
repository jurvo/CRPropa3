#ifndef CRPROPA_ADVECTIONFIELDFROMLIST_H
#define CRPROPA_ADVECTIONFIELDFROMLIST_H

#include "crpropa/advectionField/AdvectionField.h"
#include "crpropa/Common.h"

namespace crpropa {

/**
 * @class AdvectionFieldFromList
 * @brief Advection field in z-direction (away from source) with radial dependence read from a table
 */

class AdvectionFieldFromList: public AdvectionField {
private:
    std::vector<double> radius, velocity;
public:
    AdvectionFieldFromList(std::string path);
    Vector3d getField(Vector3d &pos) const;
    void loadData(std::string path);    
};

} // namespace

#endif // CRPROPA_ADVECTIONFIELDFROMLIST_H