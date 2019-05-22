#ifndef CRPROPA_DENSITYGRID_H
#define CRPROPA_DENSITYGRID_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/GridTools.h"
#include "crpropa/Geometry.h"

namespace crpropa {
	/**
	@class DensityGrid
	@brief base class for a density Grid of matter */
	class DensityGrid: public Density
	{
		public:
		ref_ptr<ScalarGrid> densityGrid; // = new ScalarGrid(Vector3d(0.), 1, 1.);
		/** @param position position to get density
		@return density at position in parts/m^3 */
		double get(const Vector3d &position) const
		{
			return densityGrid->interpolate(position);
		}
		
		/** load function for a new Grid
		@param Grid new Grid*/
		void loadGrid(ref_ptr<ScalarGrid> Grid)
		{
			densityGrid = Grid;
			return;
		}
	}; //DensityGrid
	

	
	
} // namespace crpropa

#endif // CRPROPA_DENSITYGRID_H
