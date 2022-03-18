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
    explicit ObservationLikelihood(const AgentStateObservation<AGENT> &observation) { }

    // TODO: make this into a constructor
    static std::vector<AgentStateObservation<AGENT>>
    generateObservations(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        int nTimesteps = realTrajectory.nTimesteps();
        std::vector<AgentStateObservation<AGENT>> observations;
        for (int t=0; t<nTimesteps;++t) {
            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    AGENT agent(agentId);
                    int nObserved = Random::nextBinomial(realTrajectory(t,agent), pObserveIfPresent);
                    observations.emplace_back(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
                    assert(observations.back().support().isValidSolution(realTrajectory));
                }
            }
        }
        return observations;
    }

};


#endif //ABMCMC_OBSERVATIONLIKELIHOOD_H
