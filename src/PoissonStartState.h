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
#include "ABMPriorSampler.h"

// DOMAIN should be an ABM trajectory that implements index operator over States
// and implements a static coefficients(state) function that gives the
// coefficients of a linear transform from a trajectory to a given state occupation
template<class AGENT, class DOMAIN = ModelState<AGENT>>
class PoissonStartState: public ConstrainedFactorisedDistribution<DOMAIN> {
public:
    typedef AGENT agent_type;

    std::vector<double> lambdas;
    double kappa;

    explicit PoissonStartState(std::function<double(AGENT)> agentToLambda, double kappa): kappa(kappa) {
            lambdas.reserve(AGENT::domainSize);
            for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                lambdas.push_back(agentToLambda(AGENT(agentId)));
                addPoissonFactor(agentId, agentToLambda(AGENT(agentId)));
            }
            init();
    }


    explicit PoissonStartState(std::vector<double> lambdas, double kappa): lambdas(std::move(lambdas)), kappa(kappa) {
        init();
    }

    template<class OTHERDOMAIN>
    PoissonStartState<AGENT,OTHERDOMAIN> toDomain() {
        return PoissonStartState<AGENT,OTHERDOMAIN>(lambdas, kappa);
    }


//    std::function<const ModelState<AGENT> &()> modelStateSampler() const { return _modelStateSampler; }

//    ABMPriorSampler<DOMAIN> priorSampler() const { return ABMPriorSampler<DOMAIN>(modelStateSampler()); }

//    static std::function<const ModelState<AGENT> &()> modelStateSampler(std::function<double(AGENT)> &agentToLambda) {
//        return [agentToLambda, sample = ModelState<AGENT>()]() mutable -> const ModelState<AGENT> & {
//            for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//                sample[agentId] = Random::nextPoisson(agentToLambda(AGENT(agentId)));
//            }
//            return  sample;
//        };
//    }

    const DOMAIN &operator ()() const { return nextSample(); }

    const DOMAIN &nextSample() const {
        static thread_local DOMAIN sample;
        for(int i=0; i < AGENT::domainSize; ++i)
            sample[DOMAIN::indexOf(State<agent_type>(0, agent_type(i)))] = Random::nextPoisson(lambdas[i]);
        return sample;
    }


    void init() {
        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) addPoissonFactor(agentId, lambdas[agentId]);
    }

    void addPoissonFactor(int agentId, double lambda) {
        State<AGENT> state(0,agentId);
        if(lambda == 0.0) {
            this->constraints.push_back(X(DOMAIN::indexOf(state)) == (typename DOMAIN::value_type)0);
        } else {
            double logLambda = log(lambda);
            this->addFactor(
                    [logLambda, state, kappa = this->kappa](const DOMAIN &x) {
                        return PoissonStartState::widenedUnnormalisedPoisson(logLambda, kappa, x[state]);
                    },
                    {DOMAIN::indexOf(state)}
            );
        }
    }


    // A Poisson distribution that decays exponentially below zero
    // constant factor of exp(-lambda) is removed for computational efficiency
    static std::pair<double,bool> widenedUnnormalisedPoisson(double logLambda, double kappa, int occupation) {
        //if(occupation < 0) return std::pair(0.0, false); //TODO: TEST: no negative decay!!!! //std::pair(ABM::kappa*occupation, false);               // widening
//        if(occupation < 0) return std::pair(ABM::kappa*occupation - lgamma(1-occupation), false);
        if(occupation < 0) return std::pair(kappa*occupation - exp(logLambda), false);
        return std::pair(occupation*logLambda - lgamma(occupation+1) - exp(logLambda), true);    // log of Poisson
    }

    friend std::ostream &operator <<(std::ostream &out, const PoissonStartState<AGENT,DOMAIN> &startState) {
        out << startState.lambdas << std::endl;
        out << "kappa = " << startState.kappa << std::endl;
        out << static_cast<const ConstrainedFactorisedDistribution<DOMAIN> &>(startState);
        return out;
    }

//    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
//        return this->toTrajectoryDistribution(0);
//    }


    // convert to a distribution over trajectories
//    operator FactorisedDistribution<Trajectory<AGENT>>() {
//            FactorisedDistribution<Trajectory<AGENT>> trajectoryDistribution;
//            trajectoryDistribution.logFactors.reserve(AGENT::domainSize);
//            for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
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

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
            ar >> lambdas >> kappa;
            init();
        }
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
            ar << lambdas  << kappa;
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};


#endif //ABMCMC_POISSONSTARTSTATE_H
