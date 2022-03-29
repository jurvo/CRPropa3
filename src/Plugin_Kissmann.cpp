#include <crpropa/Plugin_Kissmann.h>
#include <crpropa/Common.h>

#include <cmath>
 
using namespace crpropa;
// -----------------------------------------------------------------------------------------------------------------------

SourceSNRKissmann::SourceSNRKissmann() : SourceFeature(){
}

double SourceSNRKissmann::f_r(double r) const{
	if(r>15*kpc)
		return 0;
	if (r > 10 * kpc)
		r = 10 * kpc;
		
	return pow(r / R_earth, alpha) * std::exp(- beta * (r - R_earth) / R_earth); 
}

void SourceSNRKissmann::prepareParticle(ParticleState& particle) const {
  	Random &random = Random::instance();
	double RPos;
	while (true){
		RPos = random.rand()*R_max;
		double fTest = random.rand()*frMax;
		double fR=f_r(RPos);
		if (fTest<=fR) {
			break;
		}
	}
	double phi = random.rand()*2*M_PI;
	Vector3d pos(cos(phi)*RPos, sin(phi)*RPos, 0.);
	particle.setPosition(pos);
}

// -----------------------------------------------------------------------------------------------------------------------


SourceSpiralArm::SourceSpiralArm() : SourceFeature() {
}

double SourceSpiralArm::sourceDensity(Vector3d pos) const {
    double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double phi = std::atan2(pos.y, pos.x) - 0.5 * M_PI; // convert coordinate system to Steiman model

    double nSource = 0;
    for(size_t i = 0; i < 4; i++) {
        double arm;
        if (r < r3) {
            arm = std::exp(- (r3 - r) / sRinner);
        }
        else {
            arm = std::exp( - (r - r3) / sRouter);
        }

        double phiArm = std::log(r / a[i]) / alpha[i];
        double dPhi = fabs(phi - phiArm);
        while (dPhi > M_PI) { 
            // bring dPhi in intervall [-pi, +pi]
            dPhi -= 2 * M_PI;
        }
        arm *= std::exp(- pow_integer<2>(dPhi / sPhi));
        
        arm *= std::exp(- 0.5 * pow_integer<2>(pos.z / sZ));
        nSource += arm;
    }
    return nSource;
}


void SourceSpiralArm::prepareParticle(ParticleState &particle) const {
    Random &random = Random::instance();
    Vector3d pos;
    while(true) {
        double r = random.rand(rMax);
        double phi = random.rand() * 2 * M_PI;
        pos.x = r * std::cos(phi);
        pos.y = r * std::sin(phi);
        pos.z = 2 * (0.5 - random.rand()) * zMax;

        double test = random.rand();
        double crit = sourceDensity(pos);
        if(test <= crit) {
            break;
        }
    }

    particle.setPosition(pos);
}

void SourceSpiralArm::setRMax(double r) {
    rMax = r;
}

void SourceSpiralArm::setZMax(double z) {
    zMax = z;
}

double SourceSpiralArm::getRMax() {
    return rMax;
}

double SourceSpiralArm::getZMax() {
    return zMax;
}