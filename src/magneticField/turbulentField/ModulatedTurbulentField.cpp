#include "crpropa/magneticField/turbulentField/ModulatedTurbulentField.h"
#include "crpropa/Common.h"

#include <fstream>
#include <string>

using namespace crpropa;

ModulatedTurbulentField::ModulatedTurbulentField(ref_ptr<TurbulentField> field, std::string filePath, double h) {
    setTurbulentField(field);
    loadData(filePath);
    this -> h = h;
}

void ModulatedTurbulentField::loadData(std::string filePath){
    std::ifstream infile(filePath.c_str());

    if(!infile.good())
        throw std::runtime_error("ModulatedTurbulentField::loadData could not open "+ filePath);

    std::istream *in;
    std::string line;
    in = &infile;

    double r, b;
    while (std::getline(*in, line))
    {
        std::stringstream stream(line);
        if (stream.peek() == '#')
            continue;
        
        stream >> r >> b;
        radii.push_back(r * kpc);
        brms.push_back(b * gauss);
    }
    infile.close();
}

void ModulatedTurbulentField::setTurbulentField(ref_ptr<TurbulentField> field){
    brmsField = field -> getBrms();
    this -> field = field;
}

Vector3d ModulatedTurbulentField::getField(const Vector3d &pos) const {
    double scale = getBrmsAtPosition(pos) / brmsField;
    return scale * field -> getField(pos);
}

double ModulatedTurbulentField::getBrmsAtPosition(const Vector3d &pos) const {
    double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double z = pos.z;

    double Br = interpolate(r, radii, brms);
    return Br * std::exp(- fabs(z) / h);
}

// Vector3d ModulatedTurbulentField::getField(const Vector3d &pos, double z) const {
//     return getField(pos);
// }