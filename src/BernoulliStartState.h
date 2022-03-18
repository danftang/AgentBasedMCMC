//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_BERNOULLISTARTSTATE_H
#define ABMCMC_BERNOULLISTARTSTATE_H

#include "FactoredConvexDistribution.h"
#include "ABM.h"
#include "State.h"

template<class AGENT>
class BernoulliStartState: public FactoredConvexDistribution<ABM::occupation_type> {
public:
    explicit BernoulliStartState(std::function<double(AGENT)> agentToProbability) {
        constraints.reserve(AGENT::domainSize());
        factors.reserve(AGENT::domainSize());
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            double p = agentToProbability(agentId);
            double logP = log(p);
            double logNotP = log(1.0-p);
            addFactor(
                    0 <= 1*State<AGENT>(0,agentId) <= 1,
                    [logP,logNotP](ABM::occupation_type occupation) {
                        return(occupation >= 1)?logP:logNotP;
                    }
            );
        }
    }
};


#endif //ABMCMC_BERNOULLISTARTSTATE_H
