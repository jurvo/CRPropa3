#ifndef CRPROPA_FERRIERE_H
#define CRPROPA_FERRIERE_H

#include "crpropa/massDistribution/Density.h"

#include <cmath>
#include <string>

namespace crpropa {
/**
 @class Ferriere
 @brief model of the distribution of hydrogen in the Milky Way
  Here in model Ferriere 2007
  seperated in 2 regions (inner, outer). The border is for R=3 kpc in galactocentric radius.
  model is discribed in
outer: ApJ, 497, 759
inner:	arxiv:	astro-ph/0702532
*/
class Ferriere: public DensityDistribution {
public:
	/** @class _HI
	@brief atomic Hydrogen component of the Ferriere model */
	class _HI : public Density
	{
		public:
		/** Coordinate transformation for the CentralMolecularZone region. Rotation arround z-axis such that X is the major axis and Y is the minor axis
		@param postion position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return position in local coordinates for the CMZ region
		*/
		Vector3d CMZTrafo(const Vector3d &position) const;
		/** Coordinate transformation for the galactic bulge disk region in galactic center. Rotation arround the x-axis, the y'-axis and the x''-axis. Difened with X along the major axis, Y along the minor axis and Z along the northern normal
		@param postion position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return position in local coordinates for the GB disk region
		*/
		Vector3d DISKTrafo(const Vector3d &position) const;
	
		_HI(int id)
		{
			setNucleusId(id);
		};
		
		/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
	 	@return density of atomic hydrogen in parts/m^3 */
		double get(const Vector3d &position) const;
		
	}; //HI
	
	/** @class _HII
	@brief ionised hydrogen component of the Ferriere model */
	class _HII : public Density
	{
		public:
		
		_HII(int id)
		{
			setNucleusId(id);
		};
		
		/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		 @return density of ionised hydrogen in parts/m^3 */
		double get(const Vector3d &position) const;
		
	}; //HII
	
	/* @class _H2
	@brief molecular Hydrogen component of the Ferriere model */
	class _H2: public Density
	{
		public:
		/** Coordinate transformation for the CentralMolecularZone region. Rotation arround z-axis such that X is the major axis and Y is the minor axis
		@param postion position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return position in local coordinates for the CMZ region
		*/
		Vector3d CMZTrafo(const Vector3d &position) const;
		/** Coordinate transformation for the galactic bulge disk region in galactic center. Rotation arround the x-axis, the y'-axis and the x''-axis. Difened with X along the major axis, Y along the minor axis and Z along the northern normal
		@param postion position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		@return position in local coordinates for the GB disk region
		*/
		Vector3d DISKTrafo(const Vector3d &position) const;
		
		_H2(int id)
		{
			setNucleusId(id);
		};
		
		/** @param position position in galactic coordinates with Earth at (-8.5kpc, 0, 0)
		 @return density of molecular hydrogen in parts/m^3 */
		double get(const Vector3d &position) const;
	}; //H2
	
	Ferriere()
	{
		HI = new _HI(nucleusId(1,1));
		HII = new _HII(nucleusId(1,1));
		H2 = new _H2(nucleusId(2,1));
	};
};

} // namespace

#endif  // CRPROPA_FERRIERE_H
