#ifndef CRPROPA_DIFFUSIONCOEFFICENT_H
#define CRPROPA_DIFFUSIONCOEFFICENT_H


#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"

namespace crpropa {

/** DIFFUSIONCOEFFICENT
 @class DiffusionCoefficent
 @brief Abstract base class for a space dependet diffusion coefficent.
 */
class DiffusionCoefficent: public Referenced {
public:
	virtual ~DiffusionCoefficent() {
	}
	virtual Vector3d getCoefficent(const Vector3d &position) const = 0;
	virtual double getDivergence(const Vector3d &position) const = 0;
};


/**
 @class DiffusionCoefficentList
 @brief Diffusion coefficent decorator implementing a superposition of fields.
 */
class DiffusionCoefficentList: public DiffusionCoefficent {
	std::vector<ref_ptr<DiffusionCoefficent> > coefficents;
public:
	void addCoefficent(ref_ptr<DiffusionCoefficent> coefficent);
	Vector3d getCoefficent(const Vector3d &position) const;
	double getDivergence(const Vector3d &position) const;
};


/**
 @class UniformDiffusionCoefficent
 @brief Spatial independent diffusion coefficent 
 */
class UniformDiffusionCoefficent: public DiffusionCoefficent {
	double value;
public:
	UniformDiffusionCoefficent(const Vector3d &value);
	Vector3d getCoefficent(const Vector3d &position) const;
	double getDivergence(const Vector3d &position) const;

	std::string getDescription() const;
};
} // namespace crpropa

#endif // CRPROPA_ADVECTIONFIELD_H
