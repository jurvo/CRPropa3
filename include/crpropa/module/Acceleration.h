#ifndef ACCELERATION_MODULE
#define ACCELERATION_MODULE

#include <crpropa/Candidate.h>
#include <crpropa/Module.h>
#include <crpropa/Vector3.h>
#include <crpropa/Units.h>

namespace crpropa
{
/** @addtogroup Acceleration
 *  @{
 */

/// @class StepLengthModifier
/// @brief Modifies the steplength of an acceleration module.
class StepLengthModifier : public Referenced
{
	public:
		/// Returns an update of the steplength
		virtual double modify(double steplength, Candidate* candidate) = 0;
};


/// @class AbstractAccelerationModule
/// @brief Core functionallity for acceleration by scattering with scatter centers
///  moving in a velocity field.
/// @details The velocity field is implicity implemented in the derived classes
///  for performance reasons. Models for the dependence of the step length of
///  the scatter process are set via modifiers.
class AbstractAccelerationModule : public Module
{
	double stepLength;
	std::vector<ref_ptr<StepLengthModifier> > modifiers;

	public:
		/// The parent's constructor need to be called on initialization!
		AbstractAccelerationModule(double _stepLength = 1. * parsec);
		// add a step length modifier to the model
		void add(StepLengthModifier *modifier);
		// update the candidate
		void process(Candidate *candidate) const;

		/// Returns the velocity vector of the scatter centers in the rest frame of the
		/// candidate.
		/// Needs to be implemented in inheriting classes.
		virtual Vector3d scatterCenterVelocity(Candidate *candidate) const = 0;

		/// Scatter the candidate with a center with given scatter center
		/// velocity into a random direction. Assumes that the
		/// candidate is ultra-relativistic (m = 0).
		void scatter(Candidate* candidate, const Vector3d& scatter_center_velocity) const;
};


/// @class SecondOrderFermi
/// @brief  Implements scattering with centers moving in isotropic directions.
///   All scatter centers have the same velocity.
class SecondOrderFermi : public AbstractAccelerationModule
{
	double scatterVelocity;
	std::vector<double> angle;
	std::vector<double> angleCDF;

	public:
		SecondOrderFermi(double _scatterVelocity=.1 * crpropa::c_light, double stepLength = 1. * crpropa::parsec, unsigned int size_of_pitchangle_table=10000);
		virtual crpropa::Vector3d scatterCenterVelocity(crpropa::Candidate *candidate) const;
};


/// @class DirectedFlowScattering
/// @brief Scattering in a directed flow of scatter centers.
/// @details Two of these region with different
///		velocities can be used to create first order Fermi scenario.
///		Thanks to Aritra Ghosh, Groningn University, for first work in 2017 on the
///		shock acceleration in CRPropa leading to this module.
class DirectedFlowScattering : public AbstractAccelerationModule
{
	private:
		crpropa::Vector3d __scatterVelocity;
	public:
		DirectedFlowScattering(crpropa::Vector3d scatterCenterVelocity, double stepLength=1. * parsec);
		virtual crpropa::Vector3d scatterCenterVelocity(crpropa::Candidate *candidate) const;
};


/// @class DirectedFlowOfScatterCenters
/// @brief In a directed flow, the step length depend on the direction of the particles as headon collisions are more likely than tail=on collisions - propagating against the flow is harder.
class DirectedFlowOfScatterCenters: public StepLengthModifier
{
	private:
		Vector3d __scatterVelocity;
	public:
	DirectedFlowOfScatterCenters(const Vector3d &scatterCenterVelocity);
	double modify(double steplength, Candidate* candidate);
};


/// @class QuasiLinearTheory
/// @brief Scales the steplength according to quasi linear theory.
/// @details Following quasi-linear theory [Schlickeiser1989], the mean free path \f$\lambda\f$ of a
///  particle with energy \f$E\f$ and charge \f$Z\f$ in a field with turbulence spectrum
///  \f$\frac{k}{k_{\min}}^{-q}\f$ is
///   \f[ \lambda = {\left(\frac{B}{\delta B}\right)}^2 {\left(R_G\;  k_{\min}\right)}^{1-q} R_G \equiv \lambda_0 {\left( \frac{E}{1 EeV}\frac{1}{Z} \right)}^{2-q} \f]
///  where \f$R_G = \frac{E}{B Z}\f$ is the gyro-radius of the
///  particles.
/// This class implements the rigidity dependent scaling factor used to modify
/// the base step length.
/// @par
/// \b [Schlickeiser1989] R. Schlickeiser, Cosmic-Ray Transport and Acceleration. II.
///  Cosmic Rays in Moving Cold Media with Application to Diffu-
///  sive Shock Wave Acceleration, The Astrophysical Journal 336
///  (1989) 264. doi:10.1086/167010.
class QuasiLinearTheory : public StepLengthModifier
{
	private:
		double __referenceEnergy;
		double __turbulenceIndex;
		double __minimumRigidity;

	public:
		QuasiLinearTheory(double referenecEnergy=1.*EeV, double turbulence_index=5./3, double minimumRigidity=0);
		double modify(double steplength, Candidate* candidate);
};


/**  @} */ // end of group Acceleration

} // namespace crpropa

#endif
