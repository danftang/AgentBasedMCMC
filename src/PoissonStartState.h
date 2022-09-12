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
class PoissonStartState: public FactorisedDistribution<ModelState<AGENT>> {
public:

    explicit PoissonStartState(std::function<double(AGENT)> agentToLambda) {
            this->logFactors.reserve(AGENT::domainSize());
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                addPoissonFactor(agentId, agentToLambda(AGENT(agentId)));
            }
    }


    explicit PoissonStartState(std::initializer_list<double> lambdas): PoissonStartState([&lambdas](AGENT agent) {
        return data(lambdas)[agent];
    }) {}


    ModelState<AGENT> nextSample() const {
        ModelState<AGENT> sample;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            sample[agentId] = Random::nextPoisson(lambda(agentId));
        }
        return sample;
    }


    // A Poisson distribution that decays exponentially below zero
    static std::pair<double,bool> widenedLogPoisson(double lambda, double logLambda, int occupation) {
        if(occupation < 0) return std::pair(ABM::kappa*occupation - lambda, false);               // widening
        return std::pair(occupation*logLambda - lambda - lgamma(occupation+1), true);    // log of Poisson
    }

    friend std::ostream &operator <<(std::ostream &out, const PoissonStartState<AGENT> &startState) {
        out << "{ ";
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            out << startState.lambda(agentId) << " ";
        }
        out << "}";
        return out;
    }

    using FactorisedDistribution<ModelState<AGENT>>::addFactor;

    void addPoissonFactor(int agentId, double lambda) {
        double logLambda = log(lambda);
        addFactor(
                SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &>(
                        [lambda,logLambda,agentId](const ModelState<AGENT> &modelState) {
                            return PoissonStartState<AGENT>::widenedLogPoisson(lambda, logLambda, modelState[agentId]); // log of Poisson
                        },
                        {agentId}
                )
        );
    }

    double lambda(int agentId) const {
        return -this->logFactors[agentId](ModelState<AGENT>::zero).first;
    }

    // convert to a distribution over trajectories
    operator FactorisedDistribution<Trajectory<AGENT>>() {
            FactorisedDistribution<Trajectory<AGENT>> trajectoryDistribution;
            trajectoryDistribution.logFactors.reserve(AGENT::domainSize());
            for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
                trajectoryDistribution.addFactor(getTrajectoryFactor(agentId));
            }
            return trajectoryDistribution;
    }

    SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &> getTrajectoryFactor(int agentId) {
        double l = lambda(agentId);
        double logLambda = log(l);
        return SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &>(
                [l,logLambda,agentId](const Trajectory<AGENT> &trajectory) {
                    return PoissonStartState<AGENT>::widenedLogPoisson(l, logLambda, trajectory[State<AGENT>(0,agentId)]); // log of Poisson
                },
                State<AGENT>(0,agentId).forwardOccupationDependencies()
        );
    }


    SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &> getModelStateFactor(int agentId) {
        double l = lambda(agentId);
        double logLambda = log(l);
        return SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &>(
                [l,logLambda,agentId](const ModelState<AGENT> &modelState) {
                    return PoissonStartState<AGENT>::widenedLogPoisson(l, logLambda, modelState[agentId]); // log of Poisson
                },
                {agentId}
        );
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            double lambda;
            ar >> lambda;
            addPoissonFactor(agentId, lambda);
        }
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            ar << lambda(agentId);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};


#endif //ABMCMC_POISSONSTARTSTATE_H
