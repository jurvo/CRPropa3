#ifndef CRPROPA_PIONDECAY_H
#define CRPROPA_PIONDECAY_H

#include "crpropa/Module.h"
#include "crpropa/Random.h"

namespace crpropa {
class PionDecay: public Module {
public:
	void process(Candidate* candidate) const;
};

} // namespace

#endif // !CRPROPA_PIONDECAY_H