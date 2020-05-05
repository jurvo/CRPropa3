#include "crpropa/ParticleState.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"

#include "HepPID/ParticleIDMethods.hh"
#include "kiss/logger.h"

#include <cstdlib>
#include <sstream>

namespace crpropa {

ParticleState::ParticleState(int id, double E, Vector3d pos, Vector3d dir): id(0), energy(0.), position(0.), direction(0.), pmass(0.), charge(0.), useTimePropagation(false)
{
	setId(id);
	setEnergy(E);
	setPosition(pos);
	setDirection(dir);

}

void ParticleState::setPosition(const Vector3d &pos) {
	position = pos;
}

const Vector3d &ParticleState::getPosition() const {
	return position;
}

void ParticleState::setDirection(const Vector3d &dir) {
	direction = dir / dir.getR();
}

const Vector3d &ParticleState::getDirection() const {
	return direction;
}

void ParticleState::setEnergy(double newEnergy) {
	energy = std::max(0., newEnergy); // prevent negative energies
}

double ParticleState::getEnergy() const {
	return energy;
}

double ParticleState::getRigidity() const {
	return fabs(energy / charge);
}

void ParticleState::setId(int newId) {
	id = newId;
	if (isNucleus(id)) {
		pmass = nuclearMass(id);
		charge = chargeNumber(id) * eplus;
		if (id < 0)
			charge *= -1; // anti-nucleus
	} else {
		if (abs(id) == 11)
			pmass = mass_electron;
		charge = HepPID::charge(id) * eplus;
	}
}

int ParticleState::getId() const {
	return id;
}

void ParticleState::setUseTimePropagation(bool use) {
	useTimePropagation = use;
}

bool ParticleState::getUseTimePropagation() const {
	return useTimePropagation;
}

double ParticleState::getMass() const {
	return pmass;
}

double ParticleState::getCharge() const {
	return charge;
}

double ParticleState::getLorentzFactor() const {
	return energy / (pmass * c_squared);
}

void ParticleState::setLorentzFactor(double lf) {
	lf = std::max(0., lf); // prevent negative Lorentz factors
	energy = lf * pmass * c_squared;
}

double ParticleState::getBeta() const {
	if(useTimePropagation ==false) {
		return 1;
	}
	double alpha = energy/(pmass*c_squared);
	if(alpha<=1){
		double energyGeV = energy/GeV;
		double restEnergy = pmass*c_squared/GeV;
		KISS_LOG_WARNING 
			<< "Particle with energy lower than the rest energy\n"
			<< "Particle " << id << ", "
			<< "E = " << energyGeV << " GeV, "
			<< "E0 = " <<restEnergy << " GeV, "
			<< "x = " << position / Mpc << " Mpc, "
			<< "p = " << direction << "\n";
		return 0;
	}
	return std::sqrt(1-1/(alpha*alpha));
}

Vector3d ParticleState::getVelocity() const {
	return direction * c_light * getBeta();
}

Vector3d ParticleState::getMomentum() const {
	return direction * (energy / c_light);
}

std::string ParticleState::getDescription() const {
	std::stringstream ss;
	ss << "Particle " << id << ", ";
	ss << "E = " << energy / EeV << " EeV, ";
	ss << "x = " << position / Mpc << " Mpc, ";
	ss << "p = " << direction;
	return ss.str();
}

} // namespace crpropa
