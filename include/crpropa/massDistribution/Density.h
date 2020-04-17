#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/ParticleID.h"

namespace crpropa {
/**
	@class Density
	@brief basic class for a single mass distribution of matter
	*/
	class Density: public Referenced
	{
		private:
			int NucleusId=nucleusId(1,1);	//default: HI
			
		public:
			/** @param id id for the type of distributed matter */
			void setNucleusId(int id)
			{
				NucleusId = id;
			}
			/** @param position position to get density
			@return density in parts/m^3 */
			virtual double getDensity(const Vector3d &position) const
			{
				return 0;
			}
			
			virtual double getNucleonDensity(const Vector3d &position) const
			{
				return massNumber(NucleusId)*getDensity(position);
			}

			/** @return nucleusid of the type of distributed matter*/
			int getNucleusId() const
			{ 
				return NucleusId;
			}
			
	};
	
	
	/**
	@class DensityDistribution
	@brief basic class for a Distribution of several types of Hydrogen.
	*/
	class DensityDistribution: public Referenced
	{
		public:
			Density* HI;
			Density* HII;
			Density* H2; 			
			
			DensityDistribution()
			{ HI = new Density(); HII = new Density(); H2 = new Density();};
			/** @param position position in galactic coordinates
	 		@return density in parts/m^3, summes all parts up */
			double getDensity(const Vector3d &position) const
			{
				return HI->getDensity(position) + HII->getDensity(position) + H2->getDensity(position);
			}
			/** @param position position in galactic coordinates
			 @return nucleon density in nucleons/m^3, all parts are summed up and H2 is weighted twice */
			double getNucleonDensity(const Vector3d &position) const
			{
				double n = HI->getNucleonDensity(position);
				n += HII->getNucleonDensity(position);
				n += H2->getNucleonDensity(position);
				return n;
			}
			/** @param position position in galactic coordinates
	 		@return density of atomic hydrogen in parts/m^3 */
			double getHIDensity(const Vector3d &position) const
			{
				return HI->getDensity(position);
			}
			/** @param position position in galactic coordinates 
			 @return density of ionised hydrogen in parts/m^3 */
			double getHIIDensity(const Vector3d &position) const
			{
				return HII->getDensity(position);
			}
			/** @param position position in galactic coordinates 
			 @return density of molecular hydrogen in parts/m^3 */
			double getH2Density(const Vector3d &position) const
			{
				return H2->getDensity(position);
			}
	};
	
	/**
	@class ConstantDensity
	@brief class for a constant density in all 3 hydrogen components
	*/
	class ConstantDensity: public DensityDistribution
	{
		public:
		/**
		@class SingleConstantDensity
		@brief basic class for a constant density in a single component
		*/
		class SingleConstantDensity: public Density
		{
			private:

				double Density = 0; /**<density number in parts/m^3>*/

			public:
				/** @return constant density in parts/m^3
				@param position position in galactic coordinates	
				*/
				double get(const Vector3d &position) const
				{ 
					return Density;
				}
				/** changes (constant) densitynumber to a new value 
				@param density new densitynumber */
				void setDensity(double density)
				{
					Density = density;
				}
				/** Constructor for single constant density
				@param id nucleus id for the type of matter
				@param density density in parts/m^3
				*/
				SingleConstantDensity (int id, double density){	
					setDensity(density);
                    setNucleusId(id);
				};
		};
		/** Constructor for constant density 
		@param nHI density for atomic hydrogen
		@param nHII density for ionised hydrogen
	 	@param nH2 density for molecular hydrogen
	 	*/
		ConstantDensity(double nHI, double nHII , double nH2)
		{	
			  HI = new SingleConstantDensity(nucleusId(1,1),nHI);
			  HII = new ConstantDensity::SingleConstantDensity(nucleusId(1,1),nHII);
			  H2 = new ConstantDensity::SingleConstantDensity(nucleusId(2,1),nH2);
		};
		/** change HI density number
		@param density new HI densitynumber */
		void setHI(double density)
		{
			HI = new SingleConstantDensity(nucleusId(1,1),density);
		}
		/** change HII density number
		@param density new HII densitynumber */
		void setHII(double density)
		{
			HII = new SingleConstantDensity(nucleusId(1,1),density);
		}
		/** change H2 density number
		@param density new H2 densitynumber */
		void setH2(double density)
		{
			H2= new SingleConstantDensity(nucleusId(2,1),density);
		} 
	};
	
}  // namespace crpropa

#endif  // CRPROPA_DENSITY_H
