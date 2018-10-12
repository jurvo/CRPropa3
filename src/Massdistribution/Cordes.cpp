#include "crpropa/Massdistribution/Cordes.h"

namespace crpropa {



double Cordes::getHIIDensity(const Vector3d &position) const {
	
	
	double n=0;
	
	double x=position.x/kpc;
	double y=position.y/kpc;
	double z=position.z/kpc;
	
	double R = sqrt(pow(x,2)+pow(y,2));	//radius in galactic disk
	
	n += 0.025*exp(-fabs(z)/1)*exp(-pow(R/20,2));	//galactocentric component
	n += 0.2*exp(-fabs(z)/0.15)*exp(-pow((R-4)/2,2));	//anular component

	// check if density is NAN
	// return 0 instead
	bool NaN = std::isnan(n);
	if(NaN == true){
		KISS_LOG_WARNING
			<< "\nDensity with 'nan' occured: \n"
			<< "postion = " << position << "\n"
			<< "density-model: Cordes 1991 \n"
			<< "density-type: HII (ionised) \n"
			<< "density is set to 0. \n";
			return 0.;
	}
		return n/cm;
}

double Cordes::getDensity(const Vector3d &position) const {

	return Cordes::getHIIDensity(position);
}

bool Cordes::getisforHI() {
	return isforHI;
}

bool Cordes::getisforHII() {
	return isforHII;
}

bool Cordes::getisforH2() {
	return isforH2;
}
}
