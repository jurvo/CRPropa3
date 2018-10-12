#include "crpropa/Massdistribution/ConstantDensity.h"

namespace crpropa{

double constantDensity::getDensity(const Vector3d *position) const {
	double n = 0;
				
	if(isforHI) 
		n += HIdensitynumber;
	if(isforHII)
		n += HIIdensitynumber;
	if(isforH2)
		n += H2densitynumber;
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforHII||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n tryed to get density although all density-types are deaktivated \n"
			<< "density-module: constantDensity\n"
			<< "returned 0 density\n"
			<< "please use constant Density with 0 \n";
	}

	return n;
}

double constantDensity::getHIDensity(const Vector3d *position) const {
		
	return HIdensitynumber;
}
	
double constantDensity::getHIIDensity(const Vector3d *position) const{
		
	return HIIdensitynumber;
}

double constantDensity::getH2Density(const Vector3d *position) const{

	return H2densitynumber;
}

bool constantDensity::getisforHI() {

	return isforHI;
}

bool constantDensity::getisforHII() {

	return isforHII;
}

bool constantDensity::getisforH2() {

	return isforH2;
}

void constantDensity::setHI(bool activate, double densitynumber) {

	isforHI = activate;
	HIdensitynumber = densitynumber;
}

void constantDensity::setHII(bool activate, double densitynumber) {

	isforHII = activate;
	HIIdensitynumber = densitynumber;
}

void constantDensity::setH2(bool activate, double densitynumber) {

	isforH2 = activate;
	H2densitynumber = densitynumber;
}


}//namespace 
