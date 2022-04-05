// Base class for ABM start state distributions.
// Create an instance by creating a derived class.
//
// Created by daniel on 01/04/2022.
//

#ifndef ABMCMC_STARTSTATEDISTRIBUTION_H
#define ABMCMC_STARTSTATEDISTRIBUTION_H

#include <functional>
#include "ABM.h"
#include "FactoredConvexDistribution.h"
#include "ModelState.h"

template<class AGENT>
class StartStateDistribution: public FactoredConvexDistribution<ABM::occupation_type> {
public:
    std::function<ModelState<AGENT>(const FactoredConvexDistribution<ABM::occupation_type> &)>  sampler;

    StartStateDistribution()=default;

    explicit StartStateDistribution(std::function<ModelState<AGENT>(const FactoredConvexDistribution<ABM::occupation_type> &)> Sampler): sampler(Sampler) { }

    ModelState<AGENT> nextSample() const { return sampler(*this); }
};


#endif //ABMCMC_STARTSTATEDISTRIBUTION_H
