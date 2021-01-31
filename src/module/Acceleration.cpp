#include "crpropa/module/Acceleration.h"
#include <crpropa/Common.h>
#include <crpropa/Random.h>
#include <cmath>

namespace crpropa
{

AbstractAccelerationModule::AbstractAccelerationModule(double _stepLength) : crpropa::Module(), stepLength(_stepLength)
{
}


void AbstractAccelerationModule::add(StepLengthModifier *modifier)
{
	modifiers.push_back(modifier);
}


void AbstractAccelerationModule::scatter(crpropa::Candidate *candidate, const crpropa::Vector3d &scatter_center_velocity) const
{
	// particle momentum in lab frame
	const double E = candidate->current.getEnergy();
	const crpropa::Vector3d p = candidate->current.getMomentum();

	// transform to rest frame of scatter center (p: prime)
	const double beta = scatter_center_velocity.getR() / crpropa::c_light;
	const double gamma = 1. / sqrt(1 - beta * beta);
	const double Ep = gamma * (E - scatter_center_velocity.dot(p));
	const crpropa::Vector3d pp = (p - scatter_center_velocity* E / (crpropa::c_light * crpropa::c_light)) * gamma;

	// scatter into random direction
	const crpropa::Vector3d pp_new = crpropa::Random::instance().randVector() * pp.getR();

	// transform back
	const double E_new = gamma * (Ep + scatter_center_velocity.dot(pp_new));
	const crpropa::Vector3d p_new = (pp_new + scatter_center_velocity * Ep / (crpropa::c_light * crpropa::c_light)) * gamma;

	// update candidate properties
	candidate->current.setEnergy(E_new);
	candidate->current.setDirection(p_new / p_new.getR());
}


void AbstractAccelerationModule::process(crpropa::Candidate *candidate) const
{
	double currentStepLength = stepLength;
	for (auto m : modifiers)
	{
			currentStepLength = m->modify(currentStepLength, candidate);
	}

	double step = candidate->getCurrentStep();
	while (step > 0)
	{
			double randDistance = -1. * log(crpropa::Random::instance().rand()) * currentStepLength;

			if (step < randDistance)
			{
					candidate->limitNextStep(0.1 * currentStepLength);
					return;
			}
			scatter(candidate, scatterCenterVelocity(candidate));
			step -= randDistance;
	}
}


SecondOrderFermi::SecondOrderFermi(double _scatterVelocity, double stepLength, unsigned int size_of_pitchangle_table) : AbstractAccelerationModule(stepLength), scatterVelocity(_scatterVelocity)
{
	setDescription("SecondOrderFermi Acceleration");
	angle.resize(size_of_pitchangle_table);
	angleCDF.resize(size_of_pitchangle_table);

	// have a discretized table of beamed pitch angles
	for (unsigned int i =0; i < size_of_pitchangle_table; i++)
	{
		angle[i] = i * M_PI / (size_of_pitchangle_table-1);
		angleCDF[i] = (angle[i] +scatterVelocity / crpropa::c_light * sin(angle[i])) / M_PI;
	}
}


crpropa::Vector3d SecondOrderFermi::scatterCenterVelocity(crpropa::Candidate *candidate) const
{
	size_t idx = crpropa::closestIndex(crpropa::Random::instance().rand(), angleCDF);
	crpropa::Vector3d rv = crpropa::Random::instance().randVector();
	crpropa::Vector3d rotationAxis = candidate->current.getDirection().cross(rv);

	rv = candidate->current.getDirection().getRotated(rotationAxis, M_PI - angle[idx]);
	return rv * scatterVelocity;
}


DirectedFlowOfScatterCenters::DirectedFlowOfScatterCenters(const Vector3d &scatterCenterVelocity) : __scatterVelocity(scatterCenterVelocity)
{
}


double DirectedFlowOfScatterCenters::modify(double steplength, Candidate* candidate)
{
	double directionModifier = (-1. * __scatterVelocity.dot(candidate->current.getDirection()) + c_light) / c_light;
	return steplength / directionModifier;
}


DirectedFlowScattering::DirectedFlowScattering(crpropa::Vector3d scatterCenterVelocity, double stepLength) : __scatterVelocity(scatterCenterVelocity), AbstractAccelerationModule(stepLength)
{

	// In a directed field of scatter centers, the probability to encounter a
	// scatter center depends on the direction of the candidate.
	StepLengthModifier *mod = new DirectedFlowOfScatterCenters(__scatterVelocity);
	this->add(mod);
}


crpropa::Vector3d DirectedFlowScattering::scatterCenterVelocity(crpropa::Candidate *candidate) const
{ // does not depend on candidate here.
	return __scatterVelocity;
}


QuasiLinearTheory::QuasiLinearTheory(double referenecEnergy, double turbulence_index, double minimumRigidity) : __referenceEnergy(referenecEnergy), __turbulenceIndex(turbulence_index), __minimumRigidity(minimumRigidity)
{
}


double QuasiLinearTheory::modify(double steplength, Candidate* candidate)
{
	if (candidate->current.getRigidity() < __minimumRigidity)
	{
		return steplength * std::pow(__minimumRigidity / (__referenceEnergy / eV), 2. - __turbulenceIndex);
	}
	else
	{
		return steplength * std::pow(candidate->current.getRigidity() / (__referenceEnergy / eV), 2. - __turbulenceIndex);
	}
}










} // namespace crpropa
