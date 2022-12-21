#ifndef CRPROPA_RADIALDENSITY_H
#define CRPROPA_RADIALDENSITY_H

#include "crpropa/massDistribution/Density.h"
#include "crpropa/Common.h"

#include <cmath>
#include <sstream>
#include <string>

namespace crpropa {

/*
@class RadialDensity
@brief Density module for r-dependent densities in HI, HII and H2 component.
As model 

n(x,y,z) = min(n_max, n_0 * ( |(x,y,z) - (x_0, y_0, z_0)| / r_0)^alpha)

with n_0 as density at r_0, n_max the 
x_0, y_0, z_0 as origin
*/
class RadialDensity: public Density {
private:
	// default mode: all density types set to 0 and no activ component
	double HIdensitynumber  = 0;  /**< density for atomic hydrogen at R0 */
	double HIIdensitynumber = 0;  /**< density for ioniesd hydrogen at R0 */
	double H2densitynumber  = 0;  /**< density for molecular hydrogen at R0 */

	double HImaxdensity		= 0;  /**< maximum density for molecular hydrogen*/
	double HIImaxdensity	= 0;  /**< maximum density for ioniesd hydrogen*/
	double H2maxdensity		= 0;  /**< maximum density for molecular hydrogen*/

	bool isHI = false;  /**< If true, HI is used for sum up in getDensity */
	bool isHII = false;  /**< If true, HII is used for sum up in getDensity */
	bool isH2 = false;  /**< If true, H2 is used for sum up in getDensity */

	double R0 				= 0;  /**< the normalization radius*/
	double alpha 			= -1; /**< the power index of the r-dependency, default -1 */
	Vector3d origin;			  /**< origin of the density distribution */

public:
	/** Constructor for constant density
	 @param HI density for atomic hydrogen
	 @param HImax max density for atomic hydrogen
	 @param HII density for ionised hydrogen
	 @param HIImax max density for atomic hydrogen
	 @param H2 density for molecular hydrogen
	 @param H2max max density for atomic hydrogen
	 @param R0 normalization radius
	 @param alpha power index of the model
	 @param o origin of the mass distribution
	 */
	RadialDensity(double HI, double HImax, double HII, double HIImax, double H2, double H2max, double R0, double alpha = -1, Vector3d o = Vector3d(0.));

	/** Get density at a given position.
	 @param position 	position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns Density in parts/m^3, sum up all activated parts
	 */
	double getDensity(const Vector3d &position) const;
	/** Get HI density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns (constant) density of HI in parts/m^3
	 */
	double getHIDensity(const Vector3d &position) const;
	/** Get HII density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns (constant) density of HII in parts/m^3 
	 */
	double getHIIDensity(const Vector3d &position) const;
	/** Get H2 density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns (constant) density of H2 in parts/m^3 
	*/
	double getH2Density(const Vector3d &position) const;
	/** Get density at a given position.
	 @param position position in galactic coordinates with Earth at (-8.5 kpc, 0, 0)
	 @returns number of nucleons/m^3, sum up all activated parts and weights H2 twice 
	 */
	double getNucleonDensity(const Vector3d &position) const;

	/** Status of HI -- active or not.
	 @returns Boolean flag with activation status of HI 
	 */
	bool getIsForHI();
	/** Status of HII -- active or not.
	 @returns Boolean flag with activation status of HII
	 */
	bool getIsForHII();
	/** Status of H2 -- active or not.
	 @returns Boolean flag with activation status of H2
	 */
	bool getIsForH2();
	/** Power index of r */
	double getAlpha();
	/** The normalization R*/
	double getR0();
	/** Origin of the mass distribution*/
	Vector3d getOrigin();


	/** Change HI status and the value of the density and its maximum.
	 @param activate 		 new activation status
	 @param densityNumber	 new density [in units of 1/meter ^ 3]
	 @param maxdensitynumber new maximum density [in units of 1/meter ^ 3]
	 */
	void setHI(bool activate, double densityNumber, double maxdensitynumber);
	/** Change HI status and the value of the density.
	 @param activate 		new activation status
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHI(bool activate, double densityNumber);
	/** Change the value of the HI density and its maximum.
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 @param maxdensitynumber new maximum density [in units of 1/meter ^ 3]
	 */
	void setHI(double densityNumber, double maxdensitynumber);
	/** Change HI status and keep density unaltered.
	 @param activate 		new activation status
	 */
	void setHI(bool activate);
	/** Change HI density and keep activation status unaltered
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHI(double densityNumber);

	/** Change HII status and the value of the density and its maximum.
	 @param activate 		 new activation status
	 @param densityNumber	 new density [in units of 1/meter ^ 3]
	 @param maxdensitynumber new maximum density [in units of 1/meter ^ 3]
	 */
	void setHII(bool activate, double densityNumber, double maxdensitynumber);
	/** Change HII status and the value of the density.
	 @param activate 		new activation status
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHII(bool activate, double densityNumber);
	/** Change the value of the HII density and its maximum.
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 @param maxdensitynumber new maximum density [in units of 1/meter ^ 3]
	 */
	void setHII(double densityNumber, double maxdensitynumber);
	/** Change HII status and keep density unaltered.
	 @param activate 		new activation status
	 */
	void setHII(bool activate);
	/** Change HII density and keep activation status unaltered
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setHII(double densityNumber);

	/** Change H2 status and the value of the density and its maximum.
	 @param activate 		 new activation status
	 @param densityNumber	 new density [in units of 1/meter ^ 3]
	 @param maxdensitynumber new maximum density [in units of 1/meter ^ 3]
	 */
	void setH2(bool activate, double densityNumber, double maxdensitynumber);
	/** Change H2 status and the value of the density.
	 @param activate 		new activation status
	 @param densityNumber	new density [in units of 1/meter ^ 3]
 	 */
	void setH2(bool activate, double densityNumber);
	/** Change the value of the H2 density and its maximum.
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 @param maxdensitynumber new maximum density [in units of 1/meter ^ 3]
	 */
	void setH2(double densityNumber, double maxdensitynumber);
	/** Change H2 status and keep density unaltered.
	 @param activate 		new activation status
	 */
	void setH2(bool activate);
	/** Change H2 density and keep activation status unaltered
	 @param densityNumber	new density [in units of 1/meter ^ 3]
	 */
	void setH2(double densityNumber);

	/** Change power index */
	void setAlpha(double newAlpha);
	/** Change normalization r */
	void setR0(double newR0);
	/** Change origin of the mass distribution */
	void setOrigin(Vector3d newOrigin);

	std::string getDescription();
};

}  // namespace crpropa

#endif  // CRPROPA_RADIALDENSITY_H