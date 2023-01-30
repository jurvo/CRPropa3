#include "crpropa/module/DiffusionCoefficent.h"


namespace crpropa {

// ---- DiffusionCoefficent ----
#pragma region Static members

const double DiffusionCoefficent::D0 = 6.1e24;

#pragma endregion


#pragma region DiffusionCoefficentList
/*
void DiffusionCoefficentList::addCoefficent(ref_ptr<DiffusionCoefficent> coefficent) {
	coefficents.push_back(coefficent);
}

Vector3d DiffusionCoefficentList::getCoefficent(const Vector3d &position) const {
	Vector3d b(0.);
	for (int i = 0; i < coefficents.size(); i++)
		b += coefficents[i]->getCoefficent(position);
	return b;
}

double DiffusionCoefficentList::getDivergence(const Vector3d &position) const {
	double D=0.;
	// Work on default values for divergence or an error handling
	for (int i = 0; i < coefficents.size(); i++)
		D += coefficents[i]->getDivergence(position);
	return D;
}
*/
#pragma endregion

#pragma region UniformDiffusionCoefficent

// Constructor
UniformDiffusionCoefficent::UniformDiffusionCoefficent(double epsilon, double scale, double alpha) : epsilon(epsilon), scale(scale), alpha(alpha)  { }

double UniformDiffusionCoefficent::getCoefficent() const 
{
	return scale * D0;
}

double UniformDiffusionCoefficent::getCoefficent(const Vector3d &position) const 
{
	return getCoefficent();
}

BTensor UniformDiffusionCoefficent::getBTensor(double r) const
{
	BTensor b;
	double d = scale * D0 * pow((std::abs(r) / 4.0e9), alpha);

	b.parallel = pow(2  * d, 0.5);
	b.perpendicular = pow(2 * epsilon * d, 0.5);

	return b;
}

void UniformDiffusionCoefficent::setEpsilon(double e)
{
	if ((e > 1) or (e < 0))
		throw std::runtime_error(
				"UniformDiffusionCoefficent: epsilon not in range 0-1");
	epsilon = e;
}

void UniformDiffusionCoefficent::setAlpha(double a) 
{
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"UniformDiffusionCoefficent: alpha not in range 0-2");
	alpha = a;
}

void UniformDiffusionCoefficent::setScale(double s) 
{
	if (s < 0)
		throw std::runtime_error(
				"UniformDiffusionCoefficent: Scale error: Scale < 0");
	scale = s;
}

double UniformDiffusionCoefficent::getEpsilon() const
{
	return epsilon;
}

double UniformDiffusionCoefficent::getAlpha() const
{
	return alpha;
}

double UniformDiffusionCoefficent::getScale() const
{
	return scale;
}

std::string UniformDiffusionCoefficent::getDescription() const 
{
	std::stringstream s;
	s << "UniformDiffusionCoefficent\n";
	s << "D: " << getCoefficent()  << " m^2/s";
	return s.str();
}


#pragma endregion
/*
// pos, dir, z not needed for now, but maybe later?
void UniformDiffusionCoefficent::getBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const 
{
    double DifCoeff = scale * 6.1e24 * pow((std::abs(r) / 4.0e9), alpha);

    BTen[0] = pow( 2  * DifCoeff, 0.5); // parallel
    BTen[1] = pow(2 * epsilon * DifCoeff, 0.5); // perpendicular
    return;
}
*/

/*Vector3d UniformDiffusionCoefficent::Vector3d getGradient(const Vector3d &position) const 
{
	return 0.;
}*/



// ---- tbc ----

} // namespace crpropa
