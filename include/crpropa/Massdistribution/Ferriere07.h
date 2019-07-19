#ifndef CRPROPA_FERRIERE07_H
#define CRPROPA_FERRIERE07_H

#include "crpropa/Massdistribution/Density.h"
#include "crpropa/Vector3.h"
#include "crpropa/Units.h"

#include <math.h>

#include "kiss/logger.h"


namespace crpropa {
/**
 @class Ferriere 
 @brief distribution of hydrogen in the Milky Way 
 * Here in model Ferriere 2007
 * seperated in 2 regions (inner, outer). The border is for R=3 kpc in galactocentric radius. 
 * model is discribed in 
Außen: ApJ, 497, 759
Innen:	arxiv:	astro-ph/0702532
*/

class Ferriere: public Density {

private:

	bool isforHI = true;		// standard for all kind of distribution
	bool isforHII = true;
	bool isforH2 = true;
	double Rsun = 8500*pc;	
	

public:
	Vector3d CMZTrafo(const Vector3d &position) const; 
	Vector3d DISKTrafo(const Vector3d &position) const;

	double getDensity(const Vector3d &position) const;
	double getHIDensity(const Vector3d &position) const;
	double getHIIDensity(const Vector3d &position) const;
	double getH2Density(const Vector3d &position) const;

	void setisforHI(bool HI);
	void setisforHII(bool HII);
	void setisforH2(bool H2);
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
};

}//namespace crpropa

#endif //CRPROPA_FERRIERE07_H


