#ifndef CRPROPA_DENSITYGRID_H
#define CRPROPA_DENSITYGRID_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/GridTools.h"

namespace crpropa {
	
	class DensityGrid: public Density
	{
		public:
		ref_ptr<ScalarGrid> densityGrid; // = new ScalarGrid(Vector3d(0.), 1, 1.);
		
		double get(const Vector3d &position) const
		{
			return densityGrid->interpolate(position);
		}
		
		void loadGrid(ref_ptr<ScalarGrid> Grid)
		{
			densityGrid = Grid;
			return;
		}
	};
	
	
} // namespace crpropa

#endif // CRPROPA_DENSITYGRID_H
