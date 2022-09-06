// Represents a probability distribution over model trajectories
// based on the number of agents in each state at time t=0
// The probability is the product of Poisson distributions on each state.
//
// Created by daniel on 07/04/2022.
//

#ifndef ABMCMC_POISSONSTARTSTATE_H
#define ABMCMC_POISSONSTARTSTATE_H

#include <functional>
#include "ModelState.h"
#include "FactorisedDistribution.h"

template<class AGENT>
class PoissonStartState: public FactorisedDistribution<Trajectory<AGENT>> {
public:

    explicit PoissonStartState(std::function<double(AGENT)> agentToLambda) {
            this->logFactors.reserve(AGENT::domainSize());
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                double lambda = agentToLambda(AGENT(agentId));
                double logLambda = log(lambda);

//                int actIndices[AGENT::actDomainSize()];
//                for(int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
//                    actIndices[actId] = Event<AGENT>(0,agentId, actId).id;
//                }

                this->addFactor(
                        SparseWidenedFunction<double, Trajectory<AGENT>>(
                        [lambda,logLambda,agentId](const Trajectory<AGENT> &trajectory) {
                            int occupation = 0;
                            for(int actId = 0; actId < AGENT::actDomainSize(); ++actId) {
                                occupation += trajectory[Event<AGENT>(0,agentId,actId)];
                            }
                            return PoissonStartState<AGENT>::widenedLogPoisson(lambda, logLambda, occupation); // log of Poisson
                        },
                        State<AGENT>(0,agentId).forwardOccupationDependencies()
                        )
                );
            }
    }


    explicit PoissonStartState(std::initializer_list<double> lambdas): PoissonStartState([&lambdas](AGENT agent) {
        return data(lambdas)[agent];
    }) {}


    ModelState<AGENT> sampleStartState() {
        ModelState<AGENT> sample;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            sample[agentId] = Random::nextPoisson(-this->logFactors[agentId](0));
        }
        return sample;
    }

    Trajectory<AGENT> nextSample(int nTimesteps) {
        Trajectory<AGENT> sample(nTimesteps);
        // TODO: implement this
        return sample;
    }

    // A Poisson distribution that decays exponentially below zero
    static std::pair<double,bool> widenedLogPoisson(double lambda, double logLambda, int occupation) {
        if(occupation < 0) return std::pair(ABM::kappa*occupation - lambda, false);               // widening
        return std::pair(occupation*logLambda - lambda - lgamma(occupation+1), true);    // log of Poisson
    }

    friend std::ostream &operator <<(std::ostream &out, const PoissonStartState<AGENT> &startState) {
        Trajectory<AGENT> zeroTrajectory(1);
        out << "{ ";
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            out << startState.logFactors[agentId].exactValue(zeroTrajectory) << " ";
        }
        out << "}";
        return out;
    }
};


#endif //ABMCMC_POISSONSTARTSTATE_H
