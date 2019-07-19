#ifndef CRPROPA_CONSTANTDENSITY_H
#define CRPROPA_CONSTANTDENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <math.h>

#include "kiss/logger.h"


namespace crpropa {

/*
@class constantDensity 
@brief Density module for constant densitys in HI, HII and H2 component. 
*/
class constantDensity: public Density {

private:
	double HIdensitynumber  = 0/cm;
	double HIIdensitynumber = 0/cm;
	double H2densitynumber  = 0/cm;
	
	bool isforHI = true;
	bool isforHII = false;
	bool isforH2 = false;

public:
	
	double getDensity(const Vector3d *position) const;
	
	double getHIDensity(const Vector3d *position) const;
	double getHIIDensity(const Vector3d *position) const;
	double getH2Density(const Vector3d *position) const;
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void setHI(bool activate, double densitynumber);
	void setHII(bool activate, double densitynumber);
	void setH2(bool activate, double densitynumber);
};

} //namespace

#endif //CRPROPA_CONSTANTDENSITY_H


