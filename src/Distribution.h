//
// Created by daniel on 09/08/2021.
//

#ifndef GLPKTEST_DISTRIBUTION_H
#define GLPKTEST_DISTRIBUTION_H

#include <vector>
#include "ConvexPMF.h"

// Before we know whether we want to sample from or take the log prob of this
// distribution, then we can use the more abstract Distribution object.
template<typename DOMAIN>
class Distribution {
public:
    virtual ConvexPMF<DOMAIN>       PMF() const =0;
    virtual std::function<DOMAIN()> sampler() const =0;
    virtual int                     nDimensions() const =0;
};


#endif //GLPKTEST_DISTRIBUTION_H
