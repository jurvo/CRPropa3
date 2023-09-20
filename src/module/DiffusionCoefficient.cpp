#include "crpropa/module/DiffusionCoefficient.h"

namespace crpropa {

// ---- DiffusionCoefficient ----
#pragma region Static members
const double DiffusionCoefficient::D0 = 6.1e24;
const double DiffusionCoefficient::E0 = 4.0e9;
#pragma endregion

#pragma region DiffusionCoefficient
void DiffusionCoefficient::setEpsilon(double e)
{
	if ((e > 1) or (e < 0))
		throw std::runtime_error(
				"UniformDiffusionCoefficient: epsilon not in range 0-1");
	epsilon = e;
}

void DiffusionCoefficient::setAlpha(double a) 
{
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"UniformDiffusionCoefficient: alpha not in range 0-2");
	alpha = a;
}

void DiffusionCoefficient::setScale(double s) 
{
	if (s < 0)
		throw std::runtime_error(
				"UniformDiffusionCoefficient: Scale error: Scale < 0");
	scale = s;
}

double DiffusionCoefficient::getEpsilon() const
{
	return epsilon;
}

double DiffusionCoefficient::getAlpha() const
{
	return alpha;
}

double DiffusionCoefficient::getScale() const
{
	return scale;
}
#pragma endregion

#pragma region UniformDiffusionCoefficient
// Constructor
UniformDiffusionCoefficient::UniformDiffusionCoefficient(double epsilon, double scale, double alpha)
{ 
	this->epsilon = epsilon;
	this->scale = scale;
	this->alpha = alpha;
}

Vector3d UniformDiffusionCoefficient::getDiffusionCoefficient() const 
{
	return Vector3d(scale * D0, epsilon * scale * D0, epsilon * scale * D0);
}

Vector3d UniformDiffusionCoefficient::getDiffusionCoefficient(double r, Vector3d pos) const 
{
	return getDiffusionCoefficient();
}

Vector3d UniformDiffusionCoefficient::getDerivativeOfDiffusionCoefficient(double r, Vector3d pos) const 
{
	return Vector3d(0.);
}

Vector3d UniformDiffusionCoefficient::getBTensor(double r, Vector3d pos) const
{
	double d = scale * D0 * pow((std::abs(r) / E0), alpha);
	double k_perp = pow(2 * epsilon * d, 0.5);
	return Vector3d(pow(2  * d, 0.5), k_perp, k_perp);
}

std::string UniformDiffusionCoefficient::getDescription() const 
{
	std::stringstream s;
	s << "UniformDiffusionCoefficient\n";
	s << "D: " << getDiffusionCoefficient()  << " m^2/s\n";
	s << "scale: " << getScale() << "\n";
	s << "epsilon: " << getEpsilon() << "\n";
	s << "alpha: " << getAlpha();
	return s.str();
}
#pragma endregion

#pragma region OneDimensionalDiffusionCoefficient
// Constructor
OneDimensionalDiffusionCoefficient::OneDimensionalDiffusionCoefficient(double epsilon, double scale, double alpha, double comp, double x_sh) :
																	compression_rate(comp), shock_width(x_sh)
{ 
	this->epsilon = epsilon;
	this->scale = scale;
	this->alpha = alpha;
}

Vector3d OneDimensionalDiffusionCoefficient::getDiffusionCoefficient(double r, double x) const
{
	double d_up = 1.;
	double d_down = d_up/compression_rate;

	double a = (d_up + d_down)*0.5;
	double b = (d_up - d_down)*0.5;

	double d = pow(a - b*tanh(x/shock_width), 2) * scale * D0 * pow((std::abs(r) / E0), alpha);

	return Vector3d(d, epsilon * d, epsilon * d);
}

Vector3d OneDimensionalDiffusionCoefficient::getDiffusionCoefficient(double r, Vector3d pos) const 
{
	return getDiffusionCoefficient(r, pos.x);
}

Vector3d OneDimensionalDiffusionCoefficient::getDerivativeOfDiffusionCoefficient(double r, double x) const 
{
	double d_up = 1.;
	double d_down = d_up/compression_rate;

	double a = (d_up + d_down)*0.5;
	double b = (d_up - d_down)*0.5;

	double d = -2 * b * pow(1 / cosh(x / shock_width), 2) * (a - b*tanh(x/shock_width)) / shock_width * scale * D0 * pow((std::abs(r) / E0), alpha);

	return Vector3d(d, 0., 0.);
}

Vector3d OneDimensionalDiffusionCoefficient::getDerivativeOfDiffusionCoefficient(double r, Vector3d pos) const 
{
	return getDerivativeOfDiffusionCoefficient(r, pos.x);
}

Vector3d OneDimensionalDiffusionCoefficient::getBTensor(double r, Vector3d pos) const
{
	Vector3d d = getDiffusionCoefficient(r, pos);
	double b_par = pow(2  * d.x, 0.5);
	double b_perp = pow(2  * d.y, 0.5);
	return Vector3d(b_par, b_perp, b_perp);
}

std::string OneDimensionalDiffusionCoefficient::getDescription() const 
{
	std::stringstream s;
	s << "OneDimensionalDiffusionCoefficient\n";
	s << "scale: " << getScale() << "\n";
	s << "epsilon: " << getEpsilon() << "\n";
	s << "alpha: " << getAlpha();
	return s.str();
}

void OneDimensionalDiffusionCoefficient::setCompressionRate(double s)
{
	compression_rate = s;
}

void OneDimensionalDiffusionCoefficient::setShockwidth(double w)
{
	shock_width = w;
}

double OneDimensionalDiffusionCoefficient::getCompressionRate() const
{
	return compression_rate;
}

double OneDimensionalDiffusionCoefficient::getShockwidth() const
{
	return shock_width;
}

#pragma endregion

#pragma region Archive 
/*
// pos, dir, z not needed for now, but maybe later?
void UniformDiffusionCoefficient::getBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const 
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