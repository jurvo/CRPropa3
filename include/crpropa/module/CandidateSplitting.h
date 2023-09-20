#ifndef CRPROPA_CANDIDATEPLITTING_H
#define CRPROPA_CANDIDATEPLITTING_H

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "kiss/logger.h"


namespace crpropa {

/**
@class CandidateSplitting
@brief Candidates are split into n copies when they cross specified energy bins. Weights are set accordingly.
		In case of Diffusice Shock Acceleration, splitting number can be adapted to expected spectral index to 
		compensate for the loss of particles per magnitude in energy
*/

class CandidateSplitting: public Module {
private:
	double n_split;
	std::vector<double> Ebins;

public:

	CandidateSplitting();
	
	CandidateSplitting(int n_split, double Emin, double Emax, double n_bins);
	/** Constructor
	 @param n_split 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param n_bins		Number of energy bins 
	 */
	CandidateSplitting(int n_split, double Emin, double Emax, double n_bins, bool log);
	/** Constructor
	 @param n_split 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param n_bins		Number of energy bins 
	 @param log 		Energy bins in log
	 */

	CandidateSplitting(int SpectralIndex, double Emin, int factor);
	/** Constructor
	 @param SpectralIndex    Absolute value of expected spectral index determines splitting number 
	 @param Emin 			 Minimal energy for splitting
	 @param factor		     Determines maximal energy, Emax = Emin*10^logfactor, and Ebins
	 */

	void process(Candidate *c) const;

	void setEnergyBins(double Emin, double Emax, double n_bins, bool log);

	void setNsplit(int n);

	int getNsplit() const;

	const std::vector<double>& getEnergyBins() const;

};
/** @}*/

} // namespace crpropa
#endif // CRPROPA_PARTICLESPLITTING_H