#include "crpropa/massDistribution/Nakanishi.h"
#include "crpropa/Common.h"

namespace crpropa {

double Nakanishi::_HI::getScaleheight(const Vector3d &position) const
{
	double R = sqrt(pow_integer<2>(position.x)+pow_integer<2>(position.y));	 // radius in galactic plane
	double scaleheight = 1.06*pc*(116.3 +19.3*R/kpc+4.1*pow_integer<2>(R/kpc)-0.05*pow_integer<3>(R/kpc));
	return scaleheight;
}

double Nakanishi::_HI::getPlaneDensity(const Vector3d &position) const {
	double R = sqrt(pow_integer<2>(position.x)+pow_integer<2>(position.y));	 // radius in galactic plane
	double planedensity = 0.94/ccm*(0.6*exp(-R/(2.4*kpc))+0.24*exp(-pow_integer<2>((R-9.5*kpc)/(4.8*kpc))));
	return planedensity;
}

double Nakanishi::_HI::get(const Vector3d &position) const
{
	double planedensity = getPlaneDensity(position);
	double scaleheight = getScaleheight(position);
	double n = planedensity*pow(0.5,pow_integer<2>(position.z/scaleheight));
	return n;
}

double Nakanishi::_H2::getScaleheight(const Vector3d &position) const {
	double R = sqrt(pow_integer<2>(position.x)+ pow_integer<2>(position.y));  // radius in galactic plane
	double scaleheight = 1.06*pc*( 10.8*exp(0.28*R/kpc)+42.78);
	return scaleheight;
}

double Nakanishi::_H2::getPlaneDensity(const Vector3d &position) const {
	double R = sqrt(pow_integer<2>(position.x)+pow_integer<2>(position.y));  // radius in galactic plane
	double planedensity =0.94/ccm*(11.2*exp(-R*R/(0.874*kpc*kpc)) +0.83*exp(-pow_integer<2>((R-4*kpc)/(3.2*kpc))));
	return planedensity;
}

double Nakanishi::_H2::get(const Vector3d &position) const {
	double planedensity = getPlaneDensity(position);
	double scaleheight = getScaleheight(position);
	double n = planedensity*pow(0.5,pow_integer<2>(position.z/scaleheight));
	return n;
}

} // namespace

