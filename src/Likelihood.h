// Represents the likelihood function for a set of observations of agent state
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_LIKELIHOOD_H
#define ABMCMC_LIKELIHOOD_H

#include "ABM.h"
#include "AgentStateObservation.h"
#include "Trajectory.h"
#include "ConstrainedFactorisedDistribution.h"

template<class AGENT>
class Likelihood: public FactorisedDistribution<Trajectory<AGENT>> {
public:
    std::vector<AgentStateObservation<AGENT>> observations;

    Likelihood()=default;

    // initialise with a single observation
    explicit Likelihood(const AgentStateObservation<AGENT> &observation) {
        addObservation(observation);
    }


    explicit Likelihood(const std::vector<AgentStateObservation<AGENT>> &observations) {
        reserve(observations.size());
        for(const AgentStateObservation<AGENT> &observation: observations) {
            addObservation(observation);
        }
    }


    Likelihood(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        int nTimesteps = realTrajectory.nTimesteps();
        reserve(nTimesteps * AGENT::domainSize() * pMakeObservation);
        for (int t = 0; t < nTimesteps; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    State<AGENT> state(t, agentId);
                    int nObserved = Random::nextBinomial(realTrajectory[state], pObserveIfPresent);
                    addObservation(AgentStateObservation<AGENT>(state, nObserved, pObserveIfPresent));
                }
            }
        }
    }

    void reserve(size_t size) {
        this->logFactors.reserve(size);
        observations.reserve(size);
    }


    void addObservation(const AgentStateObservation<AGENT> &observation) {
        observations.push_back(observation);
        this->addFactor(observations.back());
    }

    friend std::ostream &operator <<(std::ostream &out, const Likelihood<AGENT> &likelihood) {
        out << "{ ";
        for(const AgentStateObservation<AGENT> &observation : likelihood.observations) {
            out << observation << " ";
        }
        out << "}";
        return out;
    }

};


#endif //ABMCMC_LIKELIHOOD_H
