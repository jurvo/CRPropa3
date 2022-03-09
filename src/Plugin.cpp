#include <crpropa/Plugin.h>
#include <crpropa/Vector3.h>
#include <crpropa/Units.h>
#include <crpropa/Common.h>

#include <string>
#include <fstream>
#include <sstream>



using namespace crpropa;

// ----------------------------------------------------------------------------

SourceRadiusFromSFR::SourceRadiusFromSFR(std::string filename) : SourceFeature() {
	loadData(filename);
}

void SourceRadiusFromSFR::loadData(std::string filename) {
	std::ifstream infile(filename.c_str());

	if(!infile.good())
		throw std::runtime_error("SourceRadiusFromSFR could not read file: " + filename);
	
	std::istream *in;
	std::string line;
	in = &infile;
	double r, rMax, s, sMax, _;

	while (std::getline(*in, line))
	{
		std::stringstream stream(line);
		if (stream.peek() == '#')
			continue;
		
		stream >> r >> s >> _;
		if (s > sMax) 
			sMax = s;
		
		if ((r>rMax) && (s>0))
			rMax = r;

		radius.push_back(r * kpc);
		sfr.push_back(s);
	}
	this -> sfrMax = sMax;
	this -> rMax = rMax * kpc;
	infile.close();
}

void SourceRadiusFromSFR::prepareParticle(ParticleState& state) const {
	Random &random = Random::instance();
	Vector3d pos = state.getPosition();
	double Rpos;
	int counter = 0;
	while (true)
	{	
		Rpos = random.rand() * rMax;
		double sPos = interpolate(Rpos, radius, sfr);
		double sTest = random.rand() * sfrMax;
		if(sTest <= sPos / 1.6){
			break;
		}
		counter++;
	}
	double phi = random.rand() * 2 * M_PI;
	pos.x = cos(phi) * Rpos;
	pos.y = sin(phi) * Rpos;
	state.setPosition(pos);
}

double SourceRadiusFromSFR::getSFRMax() {
	return sfrMax;
}

double SourceRadiusFromSFR::getRMax() {
	return rMax;
}

// ------------------------------------------------------------------------------

SourceZpositionGauss::SourceZpositionGauss(double zMax, double h) : SourceFeature(), zMax(zMax), h(h) {}

void SourceZpositionGauss::prepareParticle(ParticleState& state) const {
	Vector3d pos = state.getPosition();
	double zPos;
	Random &random = Random::instance();
	double fPos, fTest;
	while (true)
	{
		zPos = (random.rand() - 0.5) * 2 * zMax;
		fPos = exp(- zPos * zPos / h / h);
		fTest = random.rand();
		if (fTest <= fPos) {
			break;
		}
	}
	pos.z = zPos;
	state.setPosition(pos);
}

double SourceZpositionGauss::getZMax() const {
	return zMax;
}

double SourceZpositionGauss::getH() const {
	return h;
}

void SourceZpositionGauss::setZMax(double z) {
	this->zMax = z;
}

void SourceZpositionGauss::setH(double h) {
	this->h = h;
}

// -----------------------------------------------------------------------------------------------------------------------

SynchrotronSelfCompton::SynchrotronSelfCompton(ref_ptr<ModulatedTurbulentField> field, double uRad)  : Module(){
	setMagneticField(field);
	setURad(uRad);
}

void SynchrotronSelfCompton::process(Candidate *c) const {
	int id = fabs(c -> current.getId());
	
	if (id != 11)
		return; // only for electrons
	
	double E = c -> current.getEnergy();
	Vector3d pos = c -> current.getPosition();
	double dT = c -> getCurrentTimeStep(); 

	double dEdT = energyLoss(pos, E);
	double dE = dEdT * dT;
	c -> current.setEnergy(E - dE);
}

double SynchrotronSelfCompton::energyLoss(Vector3d pos, double E) const {
	double B = field -> getBrmsAtPosition(pos) / gauss;
	double beta = 8e-17 * (uRad + 6e11 * pow_integer<2>(B) / 8 / M_PI) / GeV * second;
	return beta * E * E;
}

void SynchrotronSelfCompton::setMagneticField(ref_ptr<ModulatedTurbulentField> field) {
	this -> field = field;
}

void SynchrotronSelfCompton::setURad(double u) {
	uRad = u;
}

double SynchrotronSelfCompton::getURad() const {
	return uRad;
}

std::string SynchrotronSelfCompton::getDescription() const {
	std::stringstream ss;
	ss << "continues energy loss due to Synchrotron and Inverse Compton \n";
	ss << "using a photon field with energy density uRad = " << uRad / eV * ccm << "eV / ccm \n"; 
	return ss.str();
}

// -----------------------------------------------------------------------------------------------------------------------

ObserverScaleHeight::ObserverScaleHeight(double Rmax, std::string path) : ObserverFeature(){
	this -> Rmax = Rmax;
	loadData(path);
}

void ObserverScaleHeight::loadData(std::string path) {
	std::ifstream infile(path.c_str());
	if(!infile.good())
		throw std::runtime_error("could not open file: " + path);

	double r, z;
	std::istream *in;
	std::string line;
	in = &infile;
	while(std::getline(*in, line)){
		std::stringstream stream(line);
		if(stream.peek() == '#')
			continue; // coments start with "#"

		stream >> r >> z;
		rBin.push_back(r);
		zBin.push_back(z);
	}
	infile.close();
}

DetectionState ObserverScaleHeight::checkDetection(Candidate *candidate) const {
	Vector3d pos = candidate -> current.getPosition();

	double R = std::sqrt(pos.x * pos.x + pos.y * pos.y);
	if (R > Rmax)
		return DETECTED;

	double zMax = interpolate(R, rBin, zBin);

	if (fabs(pos.z) > zMax)
		return DETECTED;

	return NOTHING;
}

// -----------------------------------------------------------------------------------------------------------------------
/*
AdvectionFieldFromList::AdvectionFieldFromList(std::string path) {
    loadData(path);
}

void AdvectionFieldFromList::loadData(std::string path) {
    std::ifstream infile(path.c_str());

    if(!infile.good())
        throw std::runtime_error("AdvectionFieldFromList could not open " + path);

    std::istream *in;
    std::string line;
    in = &infile;

    double r, v;
    while (std::getline(*in, line))
    {
        std::stringstream stream(line);
        if (stream.peek() == '#')
            continue;
        
        stream >> r >> v;
        radius.push_back(r * kpc);
        velocity.push_back(v);
    }
    infile.close();
}

Vector3d AdvectionFieldFromList::getField(Vector3d &pos) const {
    double r = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double sgnZ = pos.z / fabs(pos.z);  

    double v = interpolate(r, radius, velocity) * sgnZ;
    return Vector3d(0., 0., v);
}*/

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