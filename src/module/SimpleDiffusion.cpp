#include "crpropa/module/SimpleDiffusion.h"
using namespace crpropa;

#pragma region Constructors

// Diffusion with a magnetic field that does not have curvature
SimpleDiffusion::SimpleDiffusion(ref_ptr<MagneticField> magneticField, double tolerance,
				 double minStep, double maxStep, double epsilon) :
	minStep(0)
{
  	setMagneticField(magneticField);
  	setMaximumStep(maxStep);
  	setMinimumStep(minStep);
  	setTolerance(tolerance);
  	setEpsilon(epsilon);
  	setScale(1.);
  	setAlpha(1./3.);
	}

// Diffusion with a magnetic field that does not have curvature
SimpleDiffusion::SimpleDiffusion(ref_ptr<MagneticField> magneticField, ref_ptr<AdvectionField> advectionField, double tolerance, double minStep, double maxStep, double epsilon) :
  	minStep(0)
{
	setMagneticField(magneticField);
	setAdvectionField(advectionField);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
	setTolerance(tolerance);
	setEpsilon(epsilon);
	setScale(1.);
	setAlpha(1./3.);
  	}

#pragma endregion
	
	// TODO:
	// 1. magnetfeldlinien ohne Krümmung
	// 1.1 Integrationen fallen weg
	// 1.2 reduzierung des Diffusionstensor auf 2 Komponenten (senkrecht, parallele)
	// 		BTensor vereinfachen auf 2D
	//		Tests sollten weiterhin durchlaufen

	// 2. Diffusiontensor auslagern (ähnlich zu Advektionsfeld), Kopp, 2012, Gl. 18/19
	//  	getBTensor, getDerivatveOfBTensor
	//

void SimpleDiffusion::process(Candidate *candidate) const {

    // save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;

	// h is the time step
	double h = clip(candidate->getNextStep(), minStep, maxStep) / c_light;
	Vector3d PosIn = current.getPosition();
	Vector3d DirIn = current.getDirection();

    // rectilinear propagation for neutral particles
    // If an advection field is provided the drift is also included
	if (current.getCharge() == 0) {
		Vector3d dir = current.getDirection();
		Vector3d Pos = current.getPosition();

		Vector3d LinProp(0.);
		if (advectionField){
			driftStep(Pos, LinProp, h);
		}

		current.setPosition(Pos + LinProp + dir*h*c_light);
		candidate->setCurrentStep(h * c_light);
		candidate->setNextStep(maxStep);
		return;
	}

	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();


    // Calculate the Diffusion tensor
	double BTensor[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	calculateBTensor(rig, BTensor, PosIn, DirIn, z);


    // Generate random numbers
	double eta[] = {0., 0., 0.};
	for(size_t i=0; i < 3; i++) {
	  	eta[i] =  Random::instance().randNorm();
	}

	double TStep = BTensor[0] * eta[0];
	double PStep = BTensor[4] * (std::sqrt(eta[1] * eta[1] + eta[2] * eta[2]));
	//double BStep = BTensor[8] * eta[2]; To reduce the run time, just one perpendicular vector is drawn

	Vector3d TVec(0.);
	Vector3d PVec(0.);

	Vector3d PosOut = Vector3d(0.);
	Vector3d DirOut = Vector3d(0.);

    // Normalize the tangent vector
	TVec = magneticField->getField(PosIn, z).getUnitVector(); // Attention: Redshift is not used except for the magnitude of the field.
	// What happens if magnetfic field is 0?

    // Choose a random perpendicular vector as the Normal-vector.
    // Prevent 'nan's in the NVec-vector in the case of <TVec, NVec> = 0.
	while (PVec.getR()==0.){
	  	Vector3d RandomVector = Random::instance().randVector();
	  	PVec = TVec.cross( RandomVector );
	}
	PVec = PVec.getUnitVector();

    // Calculate the Binormal-vector - Not necessary any more
	// BVec = (TVec.cross(NVec)).getUnitVector();

	// Calculate the advection step
	Vector3d LinProp(0.);
	if (advectionField){
		driftStep(PosIn, LinProp, h);
	}
	
    // Exception: If the magnetic field vanishes: Use only advection.
    // If an advection field is not provided --> rectilinear propagation.
	double tTest = TVec.getR();
	if (tTest != tTest) {
		if (advectionField){
			current.setPosition(PosIn + LinProp);
		}
		else {
			current.setPosition(PosIn + DirIn*h*c_light);
		}
	 	candidate->setCurrentStep(h*c_light);
		double newStep = 5 * h * c_light; // max new step from CK 
		newStep = clip(newStep, minStep, maxStep);
	  	candidate->setNextStep(newStep);
	  	return;
	}

    // Integration of the SDE with a Mayorama-Euler-method aka the total step
	Vector3d PO = PosIn + LinProp + (TVec * TStep + PVec * PStep) * sqrt(h);

    // Throw error message if something went wrong with propagation.
    // Deactivate candidate.
	bool NaN = std::isnan(PO.getR());
	if (NaN == true){
		  candidate->setActive(false);
		  KISS_LOG_WARNING
			<< "\nCandidate with 'nan'-position occured: \n"
		 	<< "position = " << PO << "\n"
		  	<< "PosIn = " << PosIn << "\n"
		  	<< "TVec = " << TVec << "\n"
		  	<< "TStep = " << std::abs(TStep) << "\n"
		  	<< "NVec = " << PVec << "\n"
		  	<< "NStep = " << PStep << "\n"
//		  	<< "BVec = " << BVec << "\n"
//		  	<< "BStep = " << BStep << "\n"
			<< "Candidate is deactivated!\n";
		  return;
	}

	//DirOut = (PO - PosIn - LinProp).getUnitVector(); //Advection does not change the momentum vector
	// Random direction around the tangential direction accounts for the pitch angle average.
	DirOut = Random::instance().randConeVector(TVec, M_PI/2.);
	current.setPosition(PO);
	current.setDirection(DirOut);
	candidate->setCurrentStep(h * c_light);

	double nextStep = 4 * h * c_light; // why 4? max new step from CK, but reduced to a power of two
	nextStep = clip(nextStep, minStep, maxStep); // added this line, need more research
	candidate->setNextStep(nextStep);

    	// Debugging and Testing
    	// Delete comments if additional information should be stored in candidate
	// This property "arcLength" can be interpreted as the effective arclength
	// of the propagation along a magnetic field line.

/*
	const std::string AL = "arcLength";
	if (candidate->hasProperty(AL) == false){
	  double arcLen = (TStep + NStep + BStep) * sqrt(h);
	  candidate->setProperty(AL, arcLen);
	  return;
	}
	else {
	  double arcLen = candidate->getProperty(AL);
	  arcLen += (TStep + NStep + BStep) * sqrt(h);
	  candidate->setProperty(AL, arcLen);
	}
*/

}



void SimpleDiffusion::driftStep(const Vector3d &pos, Vector3d &linProp, double h) const {
	// add here: Drift step of diff coeff
	Vector3d advField = getAdvectionFieldAtPosition(pos);
	linProp += advField * h;
	return;
}

void SimpleDiffusion::calculateBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const {

    double DifCoeff = scale * 6.1e24 * pow((std::abs(r) / 4.0e9), alpha);
    BTen[0] = pow( 2  * DifCoeff, 0.5);
    BTen[4] = pow(2 * epsilon * DifCoeff, 0.5);
    BTen[8] = pow(2 * epsilon * DifCoeff, 0.5);
    return;
}

void SimpleDiffusion::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("SimpleDiffusion: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("SimpleDiffusion: minStep > maxStep");
	minStep = min;
}

void SimpleDiffusion::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("SimpleDiffusion: maxStep < minStep");
	maxStep = max;
}


void SimpleDiffusion::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"SimpleDiffusion: tolerance error not in range 0-1");
	tolerance = tol;
}

void SimpleDiffusion::setEpsilon(double e) {
	if ((e > 1) or (e < 0))
		throw std::runtime_error(
				"SimpleDiffusion: epsilon not in range 0-1");
	epsilon = e;
}


void SimpleDiffusion::setAlpha(double a) {
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"SimpleDiffusion: alpha not in range 0-2");
	alpha = a;
}

void SimpleDiffusion::setScale(double s) {
	if (s < 0)
		throw std::runtime_error(
				"SimpleDiffusion: Scale error: Scale < 0");
	scale = s;
}

void SimpleDiffusion::setMagneticField(ref_ptr<MagneticField> f) {
	magneticField = f;
}

void SimpleDiffusion::setAdvectionField(ref_ptr<AdvectionField> f) {
	advectionField = f;
}

double SimpleDiffusion::getMinimumStep() const {
	return minStep;
}

double SimpleDiffusion::getMaximumStep() const {
	return maxStep;
}

double SimpleDiffusion::getTolerance() const {
	return tolerance;
}

double SimpleDiffusion::getEpsilon() const {
	return epsilon;
}

double SimpleDiffusion::getAlpha() const {
	return alpha;
}

double SimpleDiffusion::getScale() const {
	return scale;
}

ref_ptr<MagneticField> SimpleDiffusion::getMagneticField() const {
	return magneticField;
}

Vector3d SimpleDiffusion::getMagneticFieldAtPosition(Vector3d pos, double z) const {
	Vector3d B(0, 0, 0);
	try {
		// check if field is valid and use the field vector at the
		// position pos with the redshift z
		if (magneticField.valid())
			B = magneticField->getField(pos, z);
	}
	catch (std::exception &e) {
		KISS_LOG_ERROR 	<< "SimpleDiffusion: Exception in SimpleDiffusion::getMagneticFieldAtPosition.\n"
				<< e.what();
	}	
	return B;
}

ref_ptr<AdvectionField> SimpleDiffusion::getAdvectionField() const {
	return advectionField;
}

Vector3d SimpleDiffusion::getAdvectionFieldAtPosition(Vector3d pos) const {
	Vector3d AdvField(0.);
	try {
		// check if field is valid and use the field vector at the
		// position pos
		if (advectionField.valid())
			AdvField = advectionField->getField(pos);
	}
	catch (std::exception &e) {
		KISS_LOG_ERROR 	<< "SimpleDiffusion: Exception in SimpleDiffusion::getAdvectionFieldAtPosition.\n"
				<< e.what();
	}
	return AdvField;
}

std::string SimpleDiffusion::getDescription() const {
	std::stringstream s;
	s << "minStep: " << minStep / kpc  << " kpc, ";
	s << "maxStep: " << maxStep / kpc  << " kpc, ";
	s << "tolerance: " << tolerance << "\n";

	if (epsilon != 0.1) {
	  s << "epsilon: " << epsilon << ", ";
	  }

	if (alpha != 1./3.) {
	  s << "alpha: " << alpha << "\n";
	  }

	if (scale != 1.) {
	  s << "D_0: " << scale*6.1e24 << " m^2/s" << "\n";
	  }

	return s.str();
}
