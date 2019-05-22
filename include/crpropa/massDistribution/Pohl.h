#ifndef CRPROPA_POHL_H
#define CRPROPA_POHL_H

#include "crpropa/Geometry.h"
#include "crpropa/massDistribution/DensityGrid.h"
//#include "crpropa/Module.h"

namespace crpropa {

/**
	@class Pohl
	@brief Grid given density of the Milky Way. Grid for H2 density based on arxiv:0712.4264 Fits-file avalibal at http://www.app.physik.uni-potsdam.de/gas.html 
	Grid for HI density based on private communication with Martin Pohl.
	*/
	class Pohl: public DensityDistribution 
	{
		public:
		class _HI: public DensityGrid
		{
			private:
			// difine values for grid
			Vector3d origin = Vector3d(-20*kpc,-20*kpc,-1.5*kpc);	
			Vector3d spacing = Vector3d(100*pc, 100*pc, 25*pc);
			ParaxialBox Box = ParaxialBox(this->origin, Vector3d(400,400,120)*spacing);
			public:
			_HI()
			{
				setNucleusId(nucleusId(1,1));
				densityGrid = new ScalarGrid(origin, 400,400,120,spacing);
				loadGridFromTxt(this->densityGrid, getDataPath("Pohl_HI_Array.txt"),1/ccm); //Fileunit parts/ccm
			};
			/** @param position position to get density in galactic coordinates with earth at (-8.5kpc, 0, 0)
			@return density in parts/m^3 */
			double get(const Vector3d &position) const
			{
				// check grid conditions to not repeat the grid
				if(Box.distance(position)<0)
				{
					//set earth to (-8.5kpc, 0, 0)
					Vector3d pos = position;
					pos.x = - pos.x;
					pos.y = - pos.y;
					return densityGrid->interpolate(pos);
				}
				else{
					return 0;
				}
			}		
		}; //_HI
		
		class _H2: public DensityGrid
		{
			private:
			//difne grid values
			Vector3d origin =Vector3d(-15*kpc,-15*kpc,-500*pc);
			Vector3d spacing = Vector3d(100*pc,100*pc,25*pc);
			ParaxialBox Box = ParaxialBox(this->origin, Vector3d(300,300,40)*spacing); //Fileunit parts/ccm
			
			public:
			_H2()
			{
				setNucleusId(nucleusId(2,1));
				densityGrid = new ScalarGrid(origin, 300,300,40,spacing);
				loadGridFromTxt(this->densityGrid, getDataPath("Pohl_H2_Array.txt"),1/ccm);
			};
			/** @param position position to get density in galactic coordinates with earth at (-8.5kpc, 0, 0)
			@return density in parts/m^3 */
			double get(const Vector3d &position) const
			{
				// check grid conditions to not repeat the grid
				if(Box.distance(position)<0)
				{
					//set earth to (-8.5kpc, 0, 0)
					Vector3d pos = position;	
					pos.x = - pos.x;
					pos.y = - pos.y;
					return densityGrid->interpolate(pos);
				}
				else{
					return 0;
				}
			 }
		}; //_H2
			
		Pohl()
		{
			HI = new _HI();
			H2 = new _H2();
		};		
	}; 
	
} // namespace crpropa

#endif //CRPROPA_POHL_H
