#ifndef CRPROPA_MODULATEDTURBULENTFIELD1
#define CRPROPA_MODULATEDTURBULENTFIELD

#include "crpropa/magneticField/turbulentField/TurbulentField.h"

#include <vector>
#include <string>

namespace crpropa{

class ModulatedTurbulentField : public MagneticField
{
private:
    ref_ptr<TurbulentField> field;
    double brmsField; // rms value of the input field. (will be normed to one)
    std::vector<double> radii; // measurement of radial position
    std::vector<double> brms; // measurement of b_rms at radius
    double h; // norm values for exponential cut
public:
    ModulatedTurbulentField(ref_ptr<TurbulentField> field, std::string filePath, double h1);

    void loadData(std::string filePath);
    void setTurbulentField(ref_ptr<TurbulentField> field);
    Vector3d getField(const Vector3d &pos) const;
    //Vector3d getField(const Vector3d pos, double z) const;
    double getBrmsAtPosition(const Vector3d &pos) const;
};

} // namespace

#endif //CRPROPA_MODULATEDTURBULENTFIELD
