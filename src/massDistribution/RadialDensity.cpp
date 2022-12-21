#include "crpropa/massDistribution/RadialDensity.h"

#include "kiss/logger.h"

#include <sstream>

namespace crpropa{

RadialDensity::RadialDensity(double HI, double HImax, double HII, double HIImax, double H2, double H2max, double R0, double alpha, Vector3d o) :
	HImaxdensity(HImax), HIImaxdensity(HIImax), H2maxdensity(H2max), R0(R0), alpha(alpha), origin(o) {
	// set all types active which are not equal 0 and change number density
	if(HI!=0)
		setHI(true, HI);
	if(HII!=0)
		setHII(true, HII);
	if(H2!=0)
		setH2(true, H2);
}
/*
RadialDensity::RadialDensity(double HI, double HImax, double HII, double HIImax, double H2, double H2max, double R0, Vector3d o) :
	HImaxdensity(HImax), HIImaxdensity(HIImax), H2maxdensity(H2max), R0(R0), alpha(-1), origin(o) { }

RadialDensity::RadialDensity(double HI, double HImax, double HII, double HIImax, double H2, double H2max, double R0, double alpha) :
	HImaxdensity(HImax), HIImaxdensity(HIImax), H2maxdensity(H2max), R0(R0), alpha(alpha), origin(Vector3d(0.)) { }

RadialDensity::RadialDensity(double HI, double HImax, double HII, double HIImax, double H2, double H2max, double R0) :
	HImaxdensity(HImax), HIImaxdensity(HIImax), H2maxdensity(H2max), R0(R0), alpha(-1), origin(Vector3d(0.)) { }
*/

#pragma region Get functions
double RadialDensity::getDensity(const Vector3d &position) const {
	double n = 0;
	double s = pow((position - origin).getR() / R0, alpha);

	if(isHI)
		n += clip(HIdensitynumber * s, 0., HImaxdensity);
	if(isHII)
		n += clip(HIIdensitynumber * s, 0., HIImaxdensity);
	if(isH2)
		n += clip(H2densitynumber * s, 0., H2maxdensity);
	return n;
}

double RadialDensity::getNucleonDensity(const Vector3d &position) const {
	double n = 0;
	double s = pow((position - origin).getR() / R0, alpha);

	if(isHI)
		n += clip(HIdensitynumber * s, 0., HImaxdensity);
	if(isHII)
		n += clip(HIIdensitynumber * s, 0., HIImaxdensity);
	if(isH2)
		n += clip(2 * H2densitynumber * s, 0., 2 * H2maxdensity);
	return n;
}

double RadialDensity::getHIDensity(const Vector3d &position) const {
	return clip(HIdensitynumber * pow((position - origin).getR() / R0, alpha), 0., HImaxdensity);
}

double RadialDensity::getHIIDensity(const Vector3d &position) const{
	return clip(HIIdensitynumber * pow((position - origin).getR() / R0, alpha), 0., HIImaxdensity);;
}

double RadialDensity::getH2Density(const Vector3d &position) const{
	return clip(H2densitynumber * pow((position - origin).getR() / R0, alpha), 0., H2maxdensity);
}

bool RadialDensity::getIsForHI() {
	return isHI;
}

bool RadialDensity::getIsForHII() {
	return isHII;
}

bool RadialDensity::getIsForH2() {
	return isH2;
}

double RadialDensity::getAlpha() {
	return alpha;
}

double RadialDensity::getR0() {
	return R0;
}

Vector3d RadialDensity::getOrigin() {
	return origin;
}

#pragma endregion

#pragma region Set functions

#pragma region HI
void RadialDensity::setHI(bool activate, double densitynumber, double maxdensitynumber) {
	isHI = activate;
	HIdensitynumber = densitynumber;
	HImaxdensity = maxdensitynumber;
}

void RadialDensity::setHI(bool activate, double densitynumber) {
	setHI(activate, densitynumber, HImaxdensity);
}

void RadialDensity::setHI(double densitynumber, double maxdensitynumber) {
	setHI(isHI, densitynumber, maxdensitynumber);
}

void RadialDensity::setHI(bool activate) {
	setHI(activate, HIdensitynumber);
}

void RadialDensity::setHI(double densitynumber) {
	setHI(isHI, densitynumber);
}
#pragma endregion

#pragma region HII
void RadialDensity::setHII(bool activate, double densitynumber, double maxdensitynumber) {
	isHII = activate;
	HIIdensitynumber = densitynumber;
	HIImaxdensity = maxdensitynumber;
}

void RadialDensity::setHII(bool activate, double densitynumber) {
	setHII(activate, densitynumber, HIImaxdensity);
}

void RadialDensity::setHII(double densitynumber, double maxdensitynumber) {
	setHII(isHII, densitynumber, maxdensitynumber);
}

void RadialDensity::setHII(bool activate) {
	setHII(activate, HIIdensitynumber);
}

void RadialDensity::setHII(double densitynumber) {
	setHII(isHII, densitynumber);
}
#pragma endregion

#pragma region H2
void RadialDensity::setH2(bool activate, double densitynumber, double maxdensitynumber) {
	isH2 = activate;
	H2densitynumber = densitynumber;
	H2maxdensity = maxdensitynumber;
}

void RadialDensity::setH2(bool activate, double densitynumber) {
	setH2(activate, densitynumber, H2maxdensity);
}

void RadialDensity::setH2(double densitynumber, double maxdensitynumber) {
	setH2(isH2, densitynumber, maxdensitynumber);
}

void RadialDensity::setH2(bool activate) {
	setHII(activate, H2densitynumber);
}

void RadialDensity::setH2(double densitynumber) {
	setH2(isH2, densitynumber);
}
#pragma endregion

void RadialDensity::setAlpha(double newAlpha) {
	alpha = newAlpha;
}

void RadialDensity::setR0(double newR0) {
	R0 = newR0;
}

void RadialDensity::setOrigin(Vector3d newOrigin) {
	origin = newOrigin;
}

#pragma endregion

std::string RadialDensity::getDescription() {
	std::stringstream s;
	s << "RadialDensity:\n";
	s<< "HI component is ";
	if(!isHI)
		s<< "not ";
	s<< "active and has a density of " << HIdensitynumber/ccm << " cm^-3";
	s<< " and a maximum of " << HImaxdensity/ccm << " cm^-3.\nHII component is ";
	if(!isHII)
		s<< "not ";
	s<< "active and has a density of " << HIIdensitynumber/ccm<<" cm^-3";
	s<< " and a maximum of " << HIImaxdensity/ccm << " cm^-3.\nH2 component is ";
	if(!isH2)
		s<<"not ";
	s<< "active and has a density of " << H2densitynumber/ccm << " cm^-3";
	s<< " and a maximum of " << H2maxdensity/ccm << " cm^-3.\n";
	s<< "Power index alpha: " << alpha << "\n";
	s<< "Normalization radius R0: " << R0/pc << " pc\n";
	s<< "Origin: " << origin;
	return s.str();
}

}  // namespace crpropa
