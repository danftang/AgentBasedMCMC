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
#include "ConstrainedFactorisedDistribution.h"

// DOMAIN should be an ABM trajectory that implements index operator over States
// and implements a static coefficients(state) function that gives the
// coefficients of a linear transform from a trajectory to a given state occupation
template<class DOMAIN, class AGENT = typename DOMAIN::agent_type>
class PoissonStartState: public ConstrainedFactorisedDistribution<DOMAIN> {
public:

    std::function<const ModelState<AGENT> &()> modelStateSampler;

    explicit PoissonStartState(std::function<double(AGENT)> agentToLambda) {
            this->factors.reserve(AGENT::domainSize());
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                addPoissonFactor(agentId, agentToLambda(AGENT(agentId)));
            }
            modelStateSampler = sampler(agentToLambda);
    }


    explicit PoissonStartState(std::initializer_list<double> lambdas): PoissonStartState([&lambdas](AGENT agent) {
        return data(lambdas)[agent];
    }) {}


    static std::function<const ModelState<AGENT> &()> sampler(std::function<double(AGENT)> &agentToLambda) {
        return [agentToLambda, sample = ModelState<AGENT>()]() mutable -> const ModelState<AGENT> & {
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                sample[agentId] = Random::nextPoisson(agentToLambda(AGENT(agentId)));
            }
            return  sample;
        };
    }


    void addPoissonFactor(int agentId, double lambda) {
        State<AGENT> state(0,agentId);
        if(lambda == 0.0) {
            auto newConstraint = this->constraints.emplace_back(DOMAIN::coefficients(state), 0);
        } else {
            double logLambda = log(lambda);
            this->addFactor(
                    SparseFunction<std::pair<double, bool>, const DOMAIN &>(
                            [logLambda, state](const DOMAIN &x) {
                                return PoissonStartState::widenedUnnormalisedPoisson(logLambda, x[state]);
                            },
                            DOMAIN::coefficients(state).indices
                    )
            );
        }
    }


    // A Poisson distribution that decays exponentially below zero
    // constant factor of exp(-lambda) is removed for computational efficiency
    static std::pair<double,bool> widenedUnnormalisedPoisson(double logLambda, int occupation) {
        if(occupation < 0) return std::pair(ABM::kappa*occupation, false);               // widening
        return std::pair(occupation*logLambda - lgamma(occupation+1), true);    // log of Poisson
    }

    friend std::ostream &operator <<(std::ostream &out, const PoissonStartState<DOMAIN> &startState) {
        out << "{ ";
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            out << startState.lambda(agentId) << " ";
        }
        out << "}";
        return out;
    }

//    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
//        return this->toTrajectoryDistribution(0);
//    }


    // convert to a distribution over trajectories
//    operator FactorisedDistribution<Trajectory<AGENT>>() {
//            FactorisedDistribution<Trajectory<AGENT>> trajectoryDistribution;
//            trajectoryDistribution.logFactors.reserve(AGENT::domainSize());
//            for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
//                trajectoryDistribution.addFactor(getTrajectoryFactor(agentId));
//            }
//            return trajectoryDistribution;
//    }



//    SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &> getTrajectoryFactor(int agentId) {
//        double l = lambda(agentId);
//        double logLambda = log(l);
//        return SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &>(
//                [l,logLambda,agentId](const Trajectory<AGENT> &trajectory) {
//                    return PoissonStartState<AGENT>::widenedLogPoisson(l, logLambda, trajectory[State<AGENT>(0,agentId)]); // log of Poisson
//                },
//                State<AGENT>(0,agentId).forwardOccupationDependencies()
//        );
//    }
//
//
//    SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &> getModelStateFactor(int agentId) {
//        double l = lambda(agentId);
//        double logLambda = log(l);
//        return SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &>(
//                [l,logLambda,agentId](const ModelState<AGENT> &modelState) {
//                    return PoissonStartState<AGENT>::widenedLogPoisson(l, logLambda, modelState[agentId]); // log of Poisson
//                },
//                {agentId}
//        );
//    }

//private:
//    friend class boost::serialization::access;
//
//    template <typename Archive>
//    void load(Archive &ar, const unsigned int version) {
//        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//            double lambda;
//            ar >> lambda;
//            addPoissonFactor(agentId, lambda);
//        }
//    }
//
//    template <typename Archive>
//    void save(Archive &ar, const unsigned int version) const {
//        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//            ar << lambda(agentId);
//        }
//    }
//
//    BOOST_SERIALIZATION_SPLIT_MEMBER();
};


#endif //ABMCMC_POISSONSTARTSTATE_H
