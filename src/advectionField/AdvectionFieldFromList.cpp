#include "crpropa/advectionField/AdvectionFieldFromList.h"

#include <fstream>
#include <string>
#include <sstream>

using namespace crpropa;

AdvectionFieldFromList::AdvectionFieldFromList(std::string path) {
    loadData(path);
}

void AdvectionFieldFromList::loadData(std::string path) {
    std::ifstream infile(path.c_str());

    if(!infile.good())
        throw std::runtime_error("AdvectionFieldFromList could not open " + path);

   std::istream *in;
    std::string line;
    in = &infile;

    double r, v;
    while (std::getline(*in, line))
    {
        std::stringstream stream(line);
        if (stream.peek() == '#')
            continue;
        
        stream >> r >> v;
        radius.push_back(r * kpc);
        velocity.push_back(v);
    }
    infile.close();
}

Vector3d AdvectionFieldFromList::getField(Vector3d &pos) const {
    double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double sgnZ = pos.z / fabs(pos.z);  

    double v = interpolate(r, radius, velocity) * sgnZ;
    return Vector3d(0., 0., v);
}