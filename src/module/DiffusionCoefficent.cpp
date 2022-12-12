#include "crpropa/module/DiffusionCoefficent.h"


namespace crpropa {

// ---- DiffusionCoefficentList ----
void DiffusionCoefficentList::addCoefficent(ref_ptr<DiffusionCoefficent> coefficent) {
	coefficents.push_back(coefficent);
}

Vector3d DiffusionCoefficentList::getCoefficent(const Vector3d &position) const {
	Vector3d b(0.);
	for (int i = 0; i < coefficents.size(); i++)
		b += coefficents[i]->getCoefficent(position);
	return b;
}

double DiffusionCoefficentList::getGradient(const Vector3d &position) const {
	double D=0.;
	// Work on default values for divergence or an error handling
	for (int i = 0; i < coefficents.size(); i++)
		D += coefficents[i]->getGradient(position);
	return D;
}

// ---- UniformDiffusionCoefficent ----
UniformDiffusionCoefficent::UniformDiffusionCoefficent(const double value) :
			value(value) {
	}

Vector3d UniformDiffusionCoefficent::getCoefficent(const Vector3d &position) const {
	return value;
	}

double UniformDiffusionCoefficent::getGradient(const Vector3d &position) const {
	return 0.;
	}

std::string UniformDiffusionCoefficent::getDescription() const {
	std::stringstream s;
	"D: " << value  << " m^2/s";
	return s.str();
}

// ---- tbc ----

} // namespace crpropa
