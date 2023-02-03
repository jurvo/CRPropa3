#include "crpropa/module/DiffusionCoefficent.h"

namespace crpropa {

// ---- DiffusionCoefficent ----
#pragma region Static members
const double DiffusionCoefficent::D0 = 6.1e24;
const double DiffusionCoefficent::E0 = 4.0e9;
#pragma endregion

#pragma region DiffusionCoefficent
void DiffusionCoefficent::setEpsilon(double e)
{
	if ((e > 1) or (e < 0))
		throw std::runtime_error(
				"UniformDiffusionCoefficent: epsilon not in range 0-1");
	epsilon = e;
}

void DiffusionCoefficent::setAlpha(double a) 
{
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"UniformDiffusionCoefficent: alpha not in range 0-2");
	alpha = a;
}

void DiffusionCoefficent::setScale(double s) 
{
	if (s < 0)
		throw std::runtime_error(
				"UniformDiffusionCoefficent: Scale error: Scale < 0");
	scale = s;
}

double DiffusionCoefficent::getEpsilon() const
{
	return epsilon;
}

double DiffusionCoefficent::getAlpha() const
{
	return alpha;
}

double DiffusionCoefficent::getScale() const
{
	return scale;
}
#pragma endregion

#pragma region UniformDiffusionCoefficent
// Constructor
UniformDiffusionCoefficent::UniformDiffusionCoefficent(double epsilon, double scale, double alpha)
{ 
	this->epsilon = epsilon;
	this->scale = scale;
	this->alpha = alpha;
}

Vector3d UniformDiffusionCoefficent::getDiffusionCoefficent() const 
{
	return Vector3d(scale * D0, epsilon * scale * D0, epsilon * scale * D0);
}

Vector3d UniformDiffusionCoefficent::getDiffusionCoefficent(double r, Vector3d pos) const 
{
	return getDiffusionCoefficent();
}

Vector3d UniformDiffusionCoefficent::getDerivativeOfDiffusionCoefficent(double r, Vector3d pos) const 
{
	return Vector3d(0.);
}

Vector3d UniformDiffusionCoefficent::getBTensor(double r, Vector3d pos) const
{
	double d = scale * D0 * pow((std::abs(r) / E0), alpha);
	double k_perp = pow(2 * epsilon * d, 0.5);
	return Vector3d(pow(2  * d, 0.5), k_perp, k_perp);
}

std::string UniformDiffusionCoefficent::getDescription() const 
{
	std::stringstream s;
	s << "UniformDiffusionCoefficent\n";
	s << "D: " << getDiffusionCoefficent()  << " m^2/s\n";
	s << "scale: " << getScale() << "\n";
	s << "epsilon: " << getEpsilon() << "\n";
	s << "alpha: " << getAlpha();
	return s.str();
}
#pragma endregion

#pragma region OneDimensionalDiffusionCoefficent
// Constructor
OneDimensionalDiffusionCoefficent::OneDimensionalDiffusionCoefficent(double epsilon, double scale, double alpha, double comp, double x_sh) :
																	compression_rate(comp), shock_width(x_sh)
{ 
	this->epsilon = epsilon;
	this->scale = scale;
	this->alpha = alpha;
}

Vector3d OneDimensionalDiffusionCoefficent::getDiffusionCoefficent(double r, double x) const
{
	double d_up = 1.;
	double d_down = d_up/compression_rate;

	double a = (d_up + d_down)*0.5;
	double b = (d_up - d_down)*0.5;

	double d = pow(a - b*tanh(x/shock_width), 2) * scale * D0 * pow((std::abs(r) / E0), alpha);

	return Vector3d(d, epsilon * d, epsilon * d);
}

Vector3d OneDimensionalDiffusionCoefficent::getDiffusionCoefficent(double r, Vector3d pos) const 
{
	return getDiffusionCoefficent(r, pos.x);
}

Vector3d OneDimensionalDiffusionCoefficent::getDerivativeOfDiffusionCoefficent(double r, double x) const 
{
	double d_up = 1.;
	double d_down = d_up/compression_rate;

	double a = (d_up + d_down)*0.5;
	double b = (d_up - d_down)*0.5;

	double d = -2 * b * pow(1 / cosh(x / shock_width), 2) * (a - b*tanh(x/shock_width)) / shock_width * scale * D0 * pow((std::abs(r) / E0), alpha);

	return Vector3d(d, 0., 0.);
}

Vector3d OneDimensionalDiffusionCoefficent::getDerivativeOfDiffusionCoefficent(double r, Vector3d pos) const 
{
	return getDerivativeOfDiffusionCoefficent(r, pos.x);
}

Vector3d OneDimensionalDiffusionCoefficent::getBTensor(double r, Vector3d pos) const
{
	Vector3d d = getDiffusionCoefficent(r, pos);
	double b_par = pow(2  * d.x, 0.5);
	double b_perp = pow(2  * d.y, 0.5);
	return Vector3d(b_par, b_perp, b_perp);
}

std::string OneDimensionalDiffusionCoefficent::getDescription() const 
{
	std::stringstream s;
	s << "OneDimensionalDiffusionCoefficent\n";
	s << "scale: " << getScale() << "\n";
	s << "epsilon: " << getEpsilon() << "\n";
	s << "alpha: " << getAlpha();
	return s.str();
}

void OneDimensionalDiffusionCoefficent::setCompressionRate(double s)
{
	compression_rate = s;
}

void OneDimensionalDiffusionCoefficent::setShockwidth(double w)
{
	shock_width = w;
}

double OneDimensionalDiffusionCoefficent::getCompressionRate() const
{
	return compression_rate;
}

double OneDimensionalDiffusionCoefficent::getShockwidth() const
{
	return shock_width;
}

#pragma endregion

#pragma region Archive 
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
#pragma endregion

// ---- tbc ----
} // namespace crpropa