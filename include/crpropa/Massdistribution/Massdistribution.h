#ifndef CRPROPA_MASSDISTRIBUTION_H
#define CRPROPA_MASSDISTRIBUTION_H


#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Massdistribution/Density.h"

#include <string>
#include <vector>

#include "kiss/logger.h"

namespace crpropa {

class Massdistribution: public Density {

ref_ptr<Density> HIDist;
ref_ptr<Density> HIIDist;
ref_ptr<Density> H2Dist;

bool isforHI=false;
bool isforHII=false;
bool isforH2=false;

bool HIisload=false;
bool HIIisload=false;
bool H2isload=false;


public:
	Massdistribution();
	double getDensity(const Vector3d &position) const;
	void add(ref_ptr<crpropa::Density> dens);	
	
	bool getisforHI();
	bool getisforHII();
	bool getisforH2();
	
	void deaktivateHI();
	void deaktivateHII();
	void deaktivateH2();

};
	
} //namespace crpropa

#endif //CRPROPA_MASSDISTRIBUTION_H


