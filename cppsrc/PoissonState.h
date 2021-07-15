//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_POISSONSTATE_H
#define GLPKTEST_POISSONSTATE_H

#include "ModelState.h"

template<typename AGENT>
class PoissonState: public ModelState<AGENT> {
public:
    int nSamples = 0;

    PoissonState(): ModelState<AGENT>() { }

    PoissonState(std::function<double (const AGENT &)> probability) {
        nSamples = 1e6;
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            AGENT agent(agentId);
            double l = probability(agent);
            if(l > 0.0) (*this)[agent] = l*nSamples - 0.5;
        }
    }

    PoissonState &operator +=(const std::vector<AGENT> &stateVector) {
        ModelState<AGENT>::operator+=(stateVector);
        ++nSamples;
        return *this;
    }


    // Addition for aggregation of states
    PoissonState &operator +=(const ModelState<AGENT> &other) {
        ModelState<AGENT>::operator+=(other);
        ++nSamples;
        return *this;
    }

    double lambda(const AGENT &agent) const {
        return ((*this)[agent]+0.5)/nSamples;
    }


    double logProb(const ModelState<AGENT> &instance) const;
    ModelState<AGENT> sample();

};


template<typename AGENT>
double PoissonState<AGENT>::logProb(const ModelState<AGENT> &instance) const {
    double logP = 0.0;
    for(auto [agent, occupation]: instance) {
        double k = fabs(occupation);
        double l = lambda(agent);
        logP += k*log(l) - l - lgamma(k+1); // log of Poisson
    }
    return logP;
}


template<typename AGENT>
ModelState<AGENT> PoissonState<AGENT>::sample() {
    ModelState<AGENT> state;
    for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
        int occupation = Random::nextPoisson(lambda(agentId));
        if(occupation > 0) state[agentId] = occupation;
    }
    return state;
}

#endif //GLPKTEST_POISSONSTATE_H
