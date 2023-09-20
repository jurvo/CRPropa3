#ifndef CRPROPA_DIFFUSIONCOEFFICIENT_H
#define CRPROPA_DIFFUSIONCOEFFICIENT_H


#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"

namespace crpropa {

#pragma region Base class DiffusionCoefficient

/** DIFFUSIONCoefficient
 @class DiffusionCoefficient
 @brief Abstract base class for a space dependet diffusion Coefficient.
 */
class DiffusionCoefficient: public Referenced {
protected:
	double epsilon;
	double scale;
	double alpha;
public:
	virtual ~DiffusionCoefficient() {}
	virtual Vector3d getDiffusionCoefficient(double r, Vector3d pos) const = 0;
	virtual Vector3d getDerivativeOfDiffusionCoefficient(double r, Vector3d pos) const = 0;
	virtual Vector3d getBTensor(double r, Vector3d pos) const = 0;
	static const double D0; // Normalization diffusion Coefficient, unit: m^2/s
	static const double E0; // Energy normalization for the diffusion Coefficient, unit: eV

	void setEpsilon(double e);
	void setAlpha(double a);
	void setScale(double s);

	double getEpsilon() const;
	double getAlpha() const;
	double getScale() const;

	virtual std::string getDescription() const = 0;
};
#pragma endregion

#pragma region UniformDiffusionCoefficient
/**
 @class UniformDiffusionCoefficient
 @brief Spatial independent diffusion Coefficient
 */
class UniformDiffusionCoefficient: public DiffusionCoefficient 
{
public:
	/**
	 * @brief Construct a new Uniform Diffusion Coefficient object
	 * 
	 * @param epsilon 
	 * @param scale 
	 * @param alpha 
	 */
	UniformDiffusionCoefficient(double epsilon, double scale, double alpha);

	Vector3d getDiffusionCoefficient() const;
	Vector3d getDiffusionCoefficient(double r, Vector3d pos) const override;

	Vector3d getDerivativeOfDiffusionCoefficient(double r, Vector3d pos) const override;

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

#pragma region OneDimensionalDiffusionCoefficient
/**
 * @class OneDimensionalDiffusionCoefficient
 * @brief DiffusionCoefficient in x-direction with shock at x = 0 and width x_sh approximated by tanh() 
		with variable compression ratio r_comp = D_up/D_down
 */
class OneDimensionalDiffusionCoefficient: public DiffusionCoefficient
{
private:
	double compression_rate; //compression ratio of shock
	double shock_width; //shock width

public:
	/**
	 * @brief Construct a new One Dimensional Diffusion Coefficient object
	 * 
	 * @param epsilon Diffusion ratio, k_perp = epsilon * k_parallel, 0 <= epsilon <= 1
	 * @param scale Scaler for the Diffuison Coefficient
	 * @param alpha Energy power index
	 * @param comp Compression ratio of the shock. It gets squared for the diffusion Coefficient profile, to presever v_adv**2 / diff across the shock!
	 * @param x_sh Shock width
	 */
	OneDimensionalDiffusionCoefficient(double epsilon, double scale, double alpha, double comp, double x_sh);

	Vector3d getDiffusionCoefficient(double r, double x) const;

	/**
	 * @brief Get the Diffusion Coefficient object
	 * 
	 * @param r Ridgidity
	 * @param pos location of the particle
	 * @return Vector3d Diffusion Coefficient
	 */ 
	Vector3d getDiffusionCoefficient(double r, Vector3d pos) const override;

	Vector3d getDerivativeOfDiffusionCoefficient(double r, double x) const;

	Vector3d getDerivativeOfDiffusionCoefficient(double r, Vector3d pos) const override;

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
