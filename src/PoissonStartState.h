//
// Created by daniel on 07/04/2022.
//

#ifndef ABMCMC_POISSONSTARTSTATE_H
#define ABMCMC_POISSONSTARTSTATE_H

#include <functional>
#include "ModelState.h"
#include "StartStateDistribution.h"

template<class AGENT>
class PoissonStartState: public StartStateDistribution<AGENT> {
public:
    static constexpr ABM::occupation_type LowerLimit = 0; //std::numeric_limits<ABM::occupation_type>::min();
    static constexpr ABM::occupation_type UpperLimit = std::numeric_limits<ABM::occupation_type>::max();

    explicit PoissonStartState(std::function<double(AGENT)> agentToLambda):
    StartStateDistribution<AGENT>(PoissonStartState<AGENT>::nextSample) {
            this->reserve(AGENT::domainSize());
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                double lambda = agentToLambda(AGENT(agentId));
                double logLambda = log(lambda);
                this->addFactor(
                        {LowerLimit <= 1*State<AGENT>(0,agentId) <= UpperLimit},
                        [lambda,logLambda](ABM::occupation_type occupation) {
                            return occupation*logLambda - lambda - lgamma(occupation+1); // log of Poisson
                        }
                );
            }
    }


    explicit PoissonStartState(std::initializer_list<double> lambdas): PoissonStartState([&lambdas](AGENT agent) {
        return data(lambdas)[agent];
    }) {}


    static ModelState<AGENT> nextSample(const FactoredConvexDistribution<ABM::occupation_type> &distribution) {
        ModelState<AGENT> sample;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            sample[agentId] = Random::nextPoisson(-distribution.factors[agentId](0));
        }
        return sample;
    }
};


#endif //ABMCMC_POISSONSTARTSTATE_H
