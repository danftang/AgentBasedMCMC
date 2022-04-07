// Represents the prior distribution over an ABM
//
// Created by daniel on 01/04/2022.
//

#ifndef ABMCMC_PRIOR_H
#define ABMCMC_PRIOR_H

#include "ABM.h"
#include "WeightedFactoredConvexDistribution.h"
#include "StartStateDistribution.h"
#include "Forecast.h"

template<class AGENT>
class Prior: public WeightedFactoredConvexDistribution<ABM::occupation_type> {
public:
    StartStateDistribution<AGENT>       startState;
    int                                 nTimesteps;

    Prior(): nTimesteps(0) { }

    Prior(int NTimesteps, StartStateDistribution<AGENT> StartState):
        WeightedFactoredConvexDistribution<ABM::occupation_type>(Forecast<AGENT>(NTimesteps) * StartState),
        nTimesteps(NTimesteps),
        startState(std::move(StartState))
    { }

    Trajectory<AGENT> nextSample() const {
        return Forecast<AGENT>::nextSample(startState, nTimesteps);
    }

};


#endif //ABMCMC_PRIOR_H
