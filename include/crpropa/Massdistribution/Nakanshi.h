#ifndef CRPROPA_NAKANSHI_H
#define CRPROPA_NAKANSHI_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <string>
#include <math.h>

#include "kiss/logger.h"

namespace crpropa {

class Nakanishi: public Density{
/*
 @class Nakanishi
 @brief Modell for HI arXiv:astro-ph/0304338
	Modell for H2 arxiv:astro-ph/0610769
*/ 

private:
	bool isforHI;
	bool isforHII;
	bool isforH2;

public:

	Nakanishi();
	double getDensity(const Vector3d &position)const;
	double getHIDensity(const Vector3d &position)const;
	double getH2Density(const Vector3d &position)const;

	double getHIScaleheight(const Vector3d &position)const;
	double getHIPlanedensity(const Vector3d &position)const;

	double getH2Scaleheight(const Vector3d &position)const;
	double getH2Planedensity(const Vector3d &position)const;


	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void setisforHI(bool HI);
	void setisforH2(bool H2);
};

} //namespace crpropa

#endif //CRPROPA_NAKANISHI_H

		

