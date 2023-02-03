#include "crpropa/module/SimpleDiffusion.h"
using namespace crpropa;

#pragma region Constructors

// Diffusion with a magnetic field that does not have curvature
SimpleDiffusion::SimpleDiffusion(ref_ptr<MagneticField> magneticField, double tolerance,
								 double minStep, double maxStep, double epsilon, double scale, double alpha) :
	minStep(0)
{
  	setMagneticField(magneticField);
  	setMaximumStep(maxStep);
  	setMinimumStep(minStep);
  	setTolerance(tolerance);
	setDiffusionCoefficent(new UniformDiffusionCoefficent(epsilon, scale, alpha));
}

// Diffusion with a magnetic field that does not have curvature
SimpleDiffusion::SimpleDiffusion(ref_ptr<MagneticField> magneticField, ref_ptr<AdvectionField> advectionField,
								 double tolerance, double minStep, double maxStep, double epsilon, double scale, double alpha) :
  	minStep(0)
{
	setMagneticField(magneticField);
	setAdvectionField(advectionField);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
	setTolerance(tolerance);
	setDiffusionCoefficent(new UniformDiffusionCoefficent(epsilon, scale, alpha));
}

SimpleDiffusion::SimpleDiffusion(ref_ptr<MagneticField> magneticField, ref_ptr<AdvectionField> advectionField,
								 double tolerance, double minStep, double maxStep, ref_ptr<DiffusionCoefficent> D) :
  	minStep(0)
{
	setMagneticField(magneticField);
	setAdvectionField(advectionField);
	setMaximumStep(maxStep);
	setMinimumStep(minStep);
	setTolerance(tolerance);
	setDiffusionCoefficent(D);
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
	// 1-D Fall:
	// 		- Magnetfeld in x-Richtung
	//		- AdvectionField mit tanh - Sophie's Git
	//		-> Diffusionskoeffizient auch mit tanh in x
	//		-> getDerivativeBTensor für DiffCoef
	//		-> v_advektion**2 / kappa_0 = const über Shock
	//		------ done so far -------------
	//		-> Advektionsschritt + (dKappa / dx, 0, 0)
	//		- AdiabaticCooling in Modulliste (mit Advektionsfeld)



/* Code example
	
	void DiffusionSDE::driftStep(const Vector3d &pos, Vector3d &linProp, double h) const {
	Vector3d advField = getAdvectionFieldAtPosition(pos);
	double DivadvField = advectionField->getDivergence(pos);
	// Sophie: calculate divkappa depending on the advection Field and use x coordinate only:
	// Future: provide diffusionCoefficient module that specifies kappa and its derivatives
	Vector3d D_DifCoeff = 0 * scale * 2 * advField * DivadvField; 
	linProp += (D_DifCoeff + advField) * h;
	return;
}

void DiffusionSDE::calculateBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const {

    double DifCoeff = scale * 6.1e24 * pow((std::abs(r) / 4.0e9), alpha); // * 6.1e24
	// Sophie: Diffusion coefficient changes with advection Field, scale now says how advection Field and Diffusion Coefficient 
	// are related (for now in 1D) + not energy dependent in Toptyghin solution:
	//Vector3d advField = getAdvectionFieldAtPosition(pos);
	DifCoeff = scale; //* advField.x * advField.x;

	BTen[0] = pow( 2  * DifCoeff, 0.5);
    BTen[4] = pow(2 * epsilon * DifCoeff, 0.5);
    BTen[8] = pow(2 * epsilon * DifCoeff, 0.5);

    return; 

}	
	*/

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

    // Calculate the Diffusion tensor, Vector3d pos, Vector3d dir, double z
	Vector3d b = DiffCoef->getBTensor(rig, PosIn);

    // Generate random numbers
	double eta[] = {0., 0., 0.};
	for(size_t i=0; i < 3; i++) {
	  	eta[i] =  Random::instance().randNorm();
	}

	double TStep = b.x * eta[0]; // parallel step
	double PStep = b.y * (std::sqrt(eta[1] * eta[1] + eta[2] * eta[2])); // perpendicular step
	// sqrt(x**2 + y**2) -> rayleigh distribution -> weibull dist with shape k = 2
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

void SimpleDiffusion::calculateBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const 
{
	
/*  double DifCoeff = scale * 6.1e24 * pow((std::abs(r) / 4.0e9), alpha);
    BTen[0] = pow( 2  * DifCoeff, 0.5);
    BTen[4] = pow(2 * epsilon * DifCoeff, 0.5);
    BTen[8] = pow(2 * epsilon * DifCoeff, 0.5);*/	
	throw std::runtime_error("SimpleDiffusion: Deprecated.");
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

void SimpleDiffusion::setDiffusionCoefficent(ref_ptr<DiffusionCoefficent> D)
{
	DiffCoef = D;
}

void SimpleDiffusion::setEpsilon(double e) {
	DiffCoef->setEpsilon(e);
}


void SimpleDiffusion::setAlpha(double a) {
	DiffCoef->setAlpha(a);
}

void SimpleDiffusion::setScale(double s) {
	DiffCoef->setScale(s);
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
	return DiffCoef->getEpsilon();
}

double SimpleDiffusion::getAlpha() const {
	return DiffCoef->getAlpha();
}

double SimpleDiffusion::getScale() const {
	return DiffCoef->getScale();
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
	s << DiffCoef->getDescription();
	return s.str();
}
