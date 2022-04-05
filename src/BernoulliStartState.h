//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_BERNOULLISTARTSTATE_H
#define ABMCMC_BERNOULLISTARTSTATE_H

#include "State.h"
#include "StartStateDistribution.h"

template<class AGENT>
class BernoulliStartState: public StartStateDistribution<AGENT> {
public:
    explicit BernoulliStartState(std::function<double(AGENT)> agentToProbability):
    StartStateDistribution<AGENT>(BernoulliStartState<AGENT>::nextSample) {
        this->reserve(AGENT::domainSize());
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            double p = agentToProbability(agentId);
            ABM::occupation_type lowerBound = p==1.0?1:0;
            ABM::occupation_type upperBound = p==0.0?0:1;
            double logP = log(p);
            double logNotP = log(1.0-p);
            this->addFactor(
                    lowerBound <= 1*State<AGENT>(0,agentId) <= upperBound,
                    [logP,logNotP](ABM::occupation_type occupation) {
                        return(occupation >= 1)?logP:logNotP;
                    }
            );
        }
    }



    explicit BernoulliStartState(std::initializer_list<double> probs): BernoulliStartState([&probs](AGENT agent) {
        return data(probs)[agent];
    }) {}


    static ModelState<AGENT> nextSample(const FactoredConvexDistribution<ABM::occupation_type> &distribution) {
        ModelState<AGENT> sample;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            if(Random::nextDouble() < exp(distribution.factors[agentId](1))) {
                sample[agentId] = 1;
            } else {
                sample[agentId] = 0;
            }
        }
        return sample;
    }
};


#endif //ABMCMC_BERNOULLISTARTSTATE_H
