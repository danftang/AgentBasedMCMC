// Represents the likelihood function for a set of observations of agent state
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_OBSERVATIONLIKELIHOOD_H
#define ABMCMC_OBSERVATIONLIKELIHOOD_H

#include "FactoredConvexDistribution.h"
#include "ABM.h"
#include "AgentStateObservation.h"
#include "Trajectory.h"

template<class AGENT>
class ObservationLikelihood: public FactoredConvexDistribution<ABM::occupation_type> {
public:

    // initialise with a single observation
    explicit ObservationLikelihood(const AgentStateObservation<AGENT> &observation) {
        addObservation(observation);
    }

    explicit ObservationLikelihood(const std::vector<AgentStateObservation<AGENT>> &observations) {
        reserve(observations.size());
        for(const AgentStateObservation<AGENT> &observation: observations) {
            addObservation(observation);
        }
    }

    void addObservation(const AgentStateObservation<AGENT> &observation) {
        addFactor(observation.constraint(), observation.toLogProbFunction());
    }

};


#endif //ABMCMC_OBSERVATIONLIKELIHOOD_H
