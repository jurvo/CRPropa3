#include "crpropa/Massdistribution/Nakanshi.h"

namespace crpropa{ 

Nakanishi::Nakanishi(){
	isforHI = true;
	isforHII = false;
	isforH2 = true;		
}



double Nakanishi::getHIScaleheight(const Vector3d &position) const{
	
	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+pow(y,2));	//Galactic Radius		
	double scaleheight = 1.06*(116.3 +19.3*R+4.1*pow(R,2)-0.05*pow(R,3));
	return scaleheight*pc;
	}

double Nakanishi::getHIPlanedensity(const Vector3d &position) const {
	
	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+pow(y,2));
	double planedensity = 0.94*(0.6*exp(-R/2.4)+0.24*exp(-pow((R-9.5)/4.8,2)));
	return planedensity/cm;
	}


double Nakanishi::getH2Scaleheight(const Vector3d &position) const {

	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+ pow(y,2));
	double Scaleheight = 1.06e-3*( 11.2*exp(0.28*R)+42.78);
	return Scaleheight*kpc;
}

double Nakanishi::getH2Planedensity(const Vector3d &position) const {

	double x = position.x/kpc;
	double y = position.y/kpc;
	double R = sqrt(pow(x,2)+pow(y,2));
	double Planedensity = 11.2*exp(-pow(R,2)/0.874) +0.83*exp(-pow((R-4)/3.2,2));
	return 0.94*Planedensity/cm;
}

double Nakanishi::getHIDensity(const Vector3d &position) const {
	
	double n = 0; //density
	double z = position.z;
	double planedensity = getHIPlanedensity(position);
	double scaleheight = getHIScaleheight(position);
	n= planedensity*pow(0.5,pow(z/scaleheight,2));
	
	// check if density is NAN
	// return 0 instead
	bool NaN = std::isnan(n);
	if(NaN == true){
		KISS_LOG_WARNING
			<< "\nDensity with 'nan' occured:\n"
			<< "postion = " << position << "\n"
			<< "density-model: Nakanishi \n"
			<< "density-type: HI (atomic)\n"
			<< "density is set to 0. \n";
			return 0;
	}
	
	return n;
}

double Nakanishi::getH2Density(const Vector3d &position) const {
	
	double n = 0; //density
	double z = position.z;
	double plane = getH2Planedensity(position);
	double scaleheight = getH2Scaleheight(position);
	n= plane*pow(0.5,pow(z/scaleheight,2));
	
	// check if density is NAN
	// return 0 instead
	bool NaN = std::isnan(n);
	if(NaN == true){
		KISS_LOG_WARNING
			<< "\nDensity with 'nan' occured:\n"
			<< "postion = " << position << "\n"
			<< "density-model: Nakanishi \n"
			<< "density-type: H2 (molecular)\n"
			<< "density is set to 0. \n";
			return 0;
	}
	
	return n;
}
	

double Nakanishi::getDensity(const Vector3d &position) const {
	double n = 0;
	if(isforHI)
		n += getHIDensity(position);
	if(isforH2)
		n += getH2Density(position);
	
	//check if any density is activ and give warning if not
	bool anyDensityActive = isforHI||isforH2;

	if(anyDensityActive == false){
		KISS_LOG_WARNING
			<< "\n tryed to get density although all density-types are deaktivated \n"
			<< "density-module: Nakanishi\n"
			<< "returned 0 density\n"
			<< "please use constant Density with 0 \n";
	}	
	
	return n;
}

bool Nakanishi::getisforHI() {
	return isforHI;
}

bool Nakanishi::getisforHII() {
	return isforHII;
}
bool Nakanishi::getisforH2() {
	return isforH2;
}

void Nakanishi::setisforHI(bool HI) {
	isforHI = HI;
}

void Nakanishi::setisforH2(bool H2) {
	isforH2 = H2;
}

} //namespace crpropa
