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

struct BTensor
{
	double parallel;
	double perpendicular;
};



/** DIFFUSIONCOEFFICENT
 @class DiffusionCoefficent
 @brief Abstract base class for a space dependet diffusion coefficent.
 */
class DiffusionCoefficent: public Referenced {
public:
	virtual ~DiffusionCoefficent() {}
	virtual double getCoefficent(const Vector3d &position) const { return 0; }
	//virtual double getDivergence(const Vector3d &position) const = 0;
	static const double D0; // Normalization diffusion coefficent, unit: m^2/s
};


/**
 @class DiffusionCoefficentList
 @brief Diffusion coefficent decorator implementing a superposition of fields.
 */
/*
class DiffusionCoefficentList: public DiffusionCoefficent {
	std::vector<ref_ptr<DiffusionCoefficent> > coefficents;
public:
	void addCoefficent(ref_ptr<DiffusionCoefficent> coefficent);
	Vector3d getCoefficent(const Vector3d &position) const;
	double getDivergence(const Vector3d &position) const;
};
*/
/**
 @class UniformDiffusionCoefficent
 @brief Spatial independent diffusion coefficent 
 */
class UniformDiffusionCoefficent: public DiffusionCoefficent {
//private:	
	double epsilon;
	double scale;
	double alpha;
public:
	UniformDiffusionCoefficent(double epsilon, double scale, double alpha);
	double getCoefficent() const;
	double getCoefficent(const Vector3d &position) const;

	/// @brief Calculates the Diffusion Tensor for a particle
	/// @param r Ridgidity of the particle
	/// @return Diffuison Tensor parallel and perpendicular component
	BTensor getBTensor(double r) const;
	//BTensor getBTensor(double r, Vector3d pos = , Vector3d dir = 0., double z = 0.) const;
	//Vector3d getGradient(const Vector3d &position) const;

	void setEpsilon(double e);
	void setAlpha(double a);
	void setScale(double s);

	double getEpsilon() const;
	double getAlpha() const;
	double getScale() const;

	std::string getDescription() const;
};
} // namespace crpropa

#endif // CRPROPA_ADVECTIONFIELD_H
