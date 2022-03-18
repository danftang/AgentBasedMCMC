// Represents the probability distribution of a trajectory given the start state,
// where the start state is implicit in the trajectory
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_TRAJECTORYFORECASTDISTRIBUTION_H
#define ABMCMC_TRAJECTORYFORECASTDISTRIBUTION_H


#include "WeightedFactoredConvexDistribution.h"
#include "TrajectoryImportance.h"

template <class AGENT>
class TrajectoryForecastDistribution: public WeightedFactoredConvexDistribution<ABM::occupation_type ,TrajectoryImportance<AGENT>> {
public:
    explicit TrajectoryForecastDistribution(int nTimesteps) { }
};


#endif //ABMCMC_TRAJECTORYFORECASTDISTRIBUTION_H
