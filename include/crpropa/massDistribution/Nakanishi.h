#ifndef CRPROPA_NAKANISHI_H
#define CRPROPA_NAKANISHI_H

#include "crpropa/massDistribution/Density.h"
#include "kiss/logger.h"

#include <cmath>
#include <sstream>

namespace crpropa {

/**
 @class Nakanishi
 @brief zylindrical symetrical model of the density distribution of the Milkyway for atomic (HI) and molecular (H2) hydrogen
 	Modell for HI arXiv:astro-ph/0304338
	Modell for H2 arxiv:astro-ph/0610769
	fit of the models given in arXiv:1607.07886
*/
class Nakanishi: public DensityDistribution {
	
	public:
	/** @class _HI
	@brief atomic hydrogen component of the Nakanishi 2003 model */
	class _HI: public Density
	{
		private:
		/** the scaleheight over the galactic plane of atomic hydrogen is fitted by polynome of degree 3
	@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	@return scaleheight at given position */
		double getScaleheight(const Vector3d &position) const;
		/** the plane density is fittet by two exponential components with e^-R and e^-(R^2)
		@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return plane density in parts/m^3 */
		double getPlaneDensity(const Vector3d &position) const;
		
		public:		
		_HI(int id)
		{
			setNucleusId(id);
		};
		
		/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return density of atomic hydrogen in parts/m^3 */
		double get(const Vector3d &position) const;
	}; //_HI 
	
	/**
	@class _H2
	@brief molecular component of the Nakanishi 2006 model */
	class _H2: public Density
	{
		private:
		/** the scaleheight over the galactic plane of molecular hydrogen is fitted by exponential function
		@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return scaleheight at given position */
		double getScaleheight(const Vector3d &position)const;
		/** the plane density is fittet by two exponential components
		@param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return plane density in parts/m^3 */
		double getPlaneDensity(const Vector3d &position)const;
		
		public:
		_H2(int id)
		{
			setNucleusId(id);
		};
		/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return density of molecular hydrogen in parts/m^3 */
		double get(const Vector3d &position) const;
	}; //_H2
	
	Nakanishi()
	{
		HI = new _HI(nucleusId(1,1));
		H2 = new _H2(nucleusId(2,1));
	};
	
};

}	//namespace

#endif // CRPROPA_NAKANISHI_H
