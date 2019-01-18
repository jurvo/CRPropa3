#ifndef CRPROPA_CORDES_H
#define CRPROPA_CORDES_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <string>

namespace crpropa {
/**
	@class Cordes
	@brief zylindrical symetrical model of the density of ionised hydrogen (HII) of the Milkyway
	Cordes et al., 1991, Nature 353,737
	*/
class Cordes: public DensityDistribution {
	
	public:
	class _HII : public Density
	{
		public:
			/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 		@return density in parts/m^3 */
			double get(const Vector3d &position) const;
			
			_HII(int id)
			{
				setNucleusId(id);
			};
		
	};
	
	Cordes()
	{
		HII = new _HII(nucleusId(1,1));
	};
	std::string getDescription();

};

}	//namespace

#endif // CRPROPA_CORDES_H
