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

template<class AGENT>
class PoissonStartState: public ConstrainedFactorisedDistribution<ModelState<AGENT>> {
public:

    explicit PoissonStartState(std::function<double(AGENT)> agentToLambda) {
            this->factors.reserve(AGENT::domainSize());
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                addPoissonFactor(agentId, agentToLambda(AGENT(agentId)));
            }
    }


    explicit PoissonStartState(std::initializer_list<double> lambdas): PoissonStartState([&lambdas](AGENT agent) {
        return data(lambdas)[agent];
    }) {}


    std::function<const ModelState<AGENT> &()> sampler() const {
        return [this, sample = ModelState<AGENT>()]() mutable -> const ModelState<AGENT> & {
            for(const auto &constraint: this->constraints) {
                assert(constraint.constant == 0);
                sample[constraint.coefficients.indices[0]] = 0;
            }
            for(const auto &factor: this->factors) {
                sample[factor.dependencies[0]] = 1; // to extract lambda from the factor
                sample[factor.dependencies[0]] = Random::nextPoisson(exp(factor(sample).first));
            }
            return  sample;
        };
    }


//    const ModelState<AGENT> &nextSample() const {
//        static thread_local ModelState<AGENT> sample;
//        for(const auto &constraint: this->constraints) {
//            sample[constraint.coefficients.indices[0]] = constraint.constant;
//        }
//        for(const auto &factor: this->logFactors) {
//            sample[factor.dependencies[0]] = Random::nextPoisson(-factor(ModelState<AGENT>::zero).first);
//        }
//        return sample;
//    }

    void addPoissonFactor(int agentId, double lambda) {
        if(lambda == 0.0) {
            this->addConstraint(1*X(agentId) == 0);
        } else {
            double logLambda = log(lambda);
            this->addFactor(
                    SparseFunction<std::pair<double, bool>, const ModelState<AGENT> &>(
                            [logLambda, agentId](const ModelState<AGENT> &modelState) {
                                return PoissonStartState<AGENT>::widenedUnnormalisedPoisson(logLambda,modelState[agentId]);
                            },
                            {agentId}
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

    friend std::ostream &operator <<(std::ostream &out, const PoissonStartState<AGENT> &startState) {
        out << "{ ";
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            out << startState.lambda(agentId) << " ";
        }
        out << "}";
        return out;
    }

    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
        return this->toTrajectoryDistribution(0);
    }


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
