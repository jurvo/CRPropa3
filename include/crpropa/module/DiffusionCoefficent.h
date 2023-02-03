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

#pragma region Base class DiffusionCoefficent

/** DIFFUSIONCOEFFICENT
 @class DiffusionCoefficent
 @brief Abstract base class for a space dependet diffusion coefficent.
 */
class DiffusionCoefficent: public Referenced {
protected:
	double epsilon;
	double scale;
	double alpha;
public:
	virtual ~DiffusionCoefficent() {}
	virtual Vector3d getDiffusionCoefficent(double r, Vector3d pos) const = 0;
	virtual Vector3d getDerivativeOfDiffusionCoefficent(double r, Vector3d pos) const = 0;
	virtual Vector3d getBTensor(double r, Vector3d pos) const = 0;
	static const double D0; // Normalization diffusion coefficent, unit: m^2/s
	static const double E0; // Energy normalization for the diffusion coefficent, unit: eV

	void setEpsilon(double e);
	void setAlpha(double a);
	void setScale(double s);

	double getEpsilon() const;
	double getAlpha() const;
	double getScale() const;

	virtual std::string getDescription() const = 0;
};
#pragma endregion

#pragma region UniformDiffusionCoefficent
/**
 @class UniformDiffusionCoefficent
 @brief Spatial independent diffusion coefficent
 */
class UniformDiffusionCoefficent: public DiffusionCoefficent 
{
public:
	/**
	 * @brief Construct a new Uniform Diffusion Coefficent object
	 * 
	 * @param epsilon 
	 * @param scale 
	 * @param alpha 
	 */
	UniformDiffusionCoefficent(double epsilon, double scale, double alpha);

	Vector3d getDiffusionCoefficent() const;
	Vector3d getDiffusionCoefficent(double r, Vector3d pos) const override;

	Vector3d getDerivativeOfDiffusionCoefficent(double r, Vector3d pos) const override;

	/**
	 * @brief Calculates the Diffusion Tensor for a particle
	 * 
	 * @param r Ridgidity of the particle 
	 * @param pos Positon of the particle
	 * @return Vector3d Diffuison Tensor parallel and perpendicular component
	 */
	Vector3d getBTensor(double r, Vector3d pos) const override;
	
	//getBTensor(double r, Vector3d pos = , Vector3d dir = 0., double z = 0.) const;

	std::string getDescription() const;
};
#pragma endregion

#pragma region OneDimensionalDiffusionCoefficent
/**
 * @class OneDimensionalDiffusionCoefficent
 * @brief DiffusionCoefficent in x-direction with shock at x = 0 and width x_sh approximated by tanh() 
		with variable compression ratio r_comp = D_up/D_down
 */
class OneDimensionalDiffusionCoefficent: public DiffusionCoefficent
{
private:
	double compression_rate; //compression ratio of shock
	double shock_width; //shock width

public:
	/**
	 * @brief Construct a new One Dimensional Diffusion Coefficent object
	 * 
	 * @param epsilon Diffusion ratio, k_perp = epsilon * k_parallel, 0 <= epsilon <= 1
	 * @param scale Scaler for the Diffuison coefficent
	 * @param alpha Energy power index
	 * @param comp Compression ratio of the shock. It gets squared for the diffusion coefficent profile, to presever v_adv**2 / diff across the shock!
	 * @param x_sh Shock width
	 */
	OneDimensionalDiffusionCoefficent(double epsilon, double scale, double alpha, double comp, double x_sh);

	Vector3d getDiffusionCoefficent(double r, double x) const;

	/**
	 * @brief Get the Diffusion Coefficent object
	 * 
	 * @param r Ridgidity
	 * @param pos location of the particle
	 * @return Vector3d Diffusion Coefficent
	 */ 
	Vector3d getDiffusionCoefficent(double r, Vector3d pos) const override;

	Vector3d getDerivativeOfDiffusionCoefficent(double r, double x) const;

	Vector3d getDerivativeOfDiffusionCoefficent(double r, Vector3d pos) const override;

	Vector3d getBTensor(double r, Vector3d pos) const override;
	
	void setCompressionRate(double s);
	void setShockwidth(double w);

	double getCompressionRate() const;
	double getShockwidth() const;

	std::string getDescription() const;
};
#pragma endregion

} // namespace crpropa

#endif // CRPROPA_ADVECTIONFIELD_H
