#ifndef CRPROPA_DENSITY_H
#define CRPROPA_DENSITY_H

#include "crpropa/Units.h"
#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/ParticleID.h"

namespace crpropa {

	class Density: public Referenced
	{
		private:
			int NucleusId=nucleusId(1,1);	//default: HI
			
		public:
			void setNucleusId(int id)
			{
				NucleusId = id;
			}
			
			virtual double get(const Vector3d &position) const
			{
				return 0;
			}
			
			int getNucleusId() const
			{ 
				return NucleusId;
			}
			
	};
	
	class DensityDistribution: public Referenced
	{
		public:
			Density* HI;
			Density* HII;
			Density* H2; 			
			
			DensityDistribution()
			{ HI = new Density(); HII = new Density(); H2 = new Density();};
			
			double getDensity(const Vector3d &position) const
			{
				return HI->get(position) + HII->get(position) + H2->get(position);
			}
			double getNucleonDensity(const Vector3d &position) const
			{
				double n = HI->get(position)*massNumber(HI->getNucleusId());
				n += HII->get(position)*massNumber(HII->getNucleusId());
				n += H2->get(position)*massNumber(H2->getNucleusId());
				return n;
			}
			
			double getHIDensity(const Vector3d &position) const
			{
				return HI->get(position);
			}
			double getHIIDensity(const Vector3d &position) const
			{
				return HII->get(position);
			}
			double getH2Density(const Vector3d &position) const
			{
				return H2->get(position);
			}
	};
	
	
	class ConstantDensity: public DensityDistribution
	{
		public:
		class SingleConstantDensity: public Density
		{
			private:

				double Density = 0;

			public:

				double get(const Vector3d &position) const
				{ 
					return Density;
				}
				void setDensity(double density)
				{
					Density = density;
				}
				
				SingleConstantDensity (int id, double density){	
					setDensity(density);
                    setNucleusId(id);
				};
		};
		
		ConstantDensity(double nHI, double nHII , double nH2)
		{	
			  HI = new SingleConstantDensity(nucleusId(1,1),nHI);
			  HII = new ConstantDensity::SingleConstantDensity(nucleusId(1,1),nHII);
			  H2 = new ConstantDensity::SingleConstantDensity(nucleusId(2,1),nH2);
		};
		
		void setHI(double density)
		{
			HI = new SingleConstantDensity(nucleusId(1,1),density);
		}
		void setHII(double density)
		{
			HII = new SingleConstantDensity(nucleusId(1,1),density);
		}
		void setH2(double density)
		{
			H2= new SingleConstantDensity(nucleusId(2,1),density);
		} 
	};
	
}  // namespace crpropa

#endif  // CRPROPA_DENSITY_H
