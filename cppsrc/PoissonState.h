//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_POISSONSTATE_H
#define GLPKTEST_POISSONSTATE_H

#include "ModelState.h"

template<typename AGENT>
class PoissonState {
public:
    int nSamples;
    ModelState<AGENT> stateCounts;

    PoissonState(): nSamples(0), stateCounts() { }

    PoissonState(std::function<double (const AGENT &)> probability) {
        nSamples = 1e6;
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            AGENT agent(agentId);
            double l = probability(agent);
            if(l > 0.0) stateCounts[agent] = l*nSamples - 0.5;
        }
    }

    PoissonState &operator +=(const std::vector<AGENT> &stateVector) {
        stateCounts += stateVector;
        ++nSamples;
        return *this;
    }


    // Addition for aggregation of states
    PoissonState &operator +=(const ModelState<AGENT> &other) {
        stateCounts += other;
        ++nSamples;
        return *this;
    }


    double lambda(int agentId) const {
        return (stateCounts[agentId]+0.5)/nSamples;
    }


    double logProb(const ModelState<AGENT> &instance) const {
        double logP = 0.0;
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            double k = fabs(instance[agentId]);
            double l = lambda(agentId);
            logP += k*log(l) - l - lgamma(k+1); // log of Poisson
        }
        return logP;
    }


    ModelState<AGENT> sample() const {
        ModelState<AGENT> state;
        for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
            int occupation = Random::nextPoisson(lambda(agentId));
            if(occupation > 0) state[agentId] = occupation;
        }
        return state;
    }

    friend std::ostream &operator <<(std::ostream &out, const PoissonState<AGENT> &poissonState) {
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            if(poissonState.stateCounts[agentId] != 0.0) out << AGENT(agentId) << " -> " << poissonState.lambda(agentId) << " ";
        }
        return  out;
    }

};

#endif //GLPKTEST_POISSONSTATE_H
