//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_BERNOULLISTARTSTATE_H
#define ABMCMC_BERNOULLISTARTSTATE_H

#include "ConstrainedFactorisedDistribution.h"
#include "ModelState.h"
#include "ABMPriorSampler.h"

template<class DOMAIN, class AGENT = typename DOMAIN::agent_type>
class BernoulliStartState: public ConstrainedFactorisedDistribution<DOMAIN> {
protected:
    std::function<const ModelState<AGENT> &()> _modelStateSampler;
public:

    using typename ConstrainedFactorisedDistribution<DOMAIN>::coefficient_type;


    explicit BernoulliStartState(std::function<double(AGENT)> agentToProb) {
        this->factors.reserve(AGENT::domainSize);
        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
            addBernoulliFactor(agentId, agentToProb(AGENT(agentId)));
        }
        _modelStateSampler = modelStateSampler(agentToProb);
    }

    explicit BernoulliStartState(std::initializer_list<double> probs): BernoulliStartState([probs = std::vector(probs)](AGENT agent) {
        return probs[agent];
    }) {}

    std::function<const ModelState<AGENT> &()> modelStateSampler() const { return _modelStateSampler; }

    ABMPriorSampler<DOMAIN> priorSampler() const { return ABMPriorSampler<DOMAIN>(modelStateSampler()); }

    static std::function<const ModelState<AGENT> &()> modelStateSampler(std::function<double(AGENT)> &agentToprob) {
        return [sample = ModelState<AGENT>(), agentToprob]() mutable -> const ModelState<AGENT> & {
            for(int agentId=0; agentId < AGENT::domainSize; ++agentId) {
                sample[agentId] = Random::nextBool(agentToprob(AGENT(agentId)));
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
//            sample[factor.dependencies[0]] = Random::nextBool(1.0-factor(ModelState<AGENT>::zero).first)?0:1;
//        }
//        return sample;
//    }


    void addBernoulliFactor(int agentId, double p) {
        State<AGENT> state(0,agentId);
        if(p != 1.0 && p != 0.0) {
            double logP = log(p);
            double logNotP = log(1.0-p);
            this->addFactor(
                    [logP, logNotP, state](const DOMAIN &x) {
                        return widenedLogBernoulli(logP, logNotP, x[state]); // log of Bernoulli
                    },
                    { DOMAIN::indexOf(state) }

            );
        } else {
            this->constraints.push_back(X(DOMAIN::indexOf(state)) == (typename DOMAIN::value_type)p);
        }
    }

//    // convert to a distribution over trajectories
//    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
//        return this->toTrajectoryDistribution(0);
//    }


//    SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &> getTrajectoryFactor(int agentId) {
//        double p = probability(agentId);
//        return SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &>(
//                [p,agentId](const Trajectory<AGENT> &trajectory) {
//                    return BernoulliStartState<AGENT>::widenedLogBernoulli(p, trajectory[State<AGENT>(0,agentId)]);
//                },
//                State<AGENT>(0,agentId).forwardOccupationDependencies()
//        );
//    }
//
//
//    SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &> getModelStateFactor(int agentId) {
//        double p = probability(agentId);
//        return SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &>(
//                [p,agentId](const ModelState<AGENT> &modelState) {
//                    return BernoulliStartState<AGENT>::widenedLogBernoulli(p, modelState[agentId]); // log of Poisson
//                },
//                {agentId}
//        );
//    }

    // A Bernoulli distribution that decays exponentially below zero
    static std::pair<double,bool> widenedLogBernoulli(double logP, double logNotP, int occupation) {
        if(occupation == 0) return std::pair(logNotP, true);
        if(occupation < 0) return std::pair(logNotP + ABM::kappa*occupation, false);
        if(occupation > 1) return std::pair(logP + ABM::kappa*(1-occupation), false);
        return std::pair(logP, true);
    }

    friend std::ostream &operator <<(std::ostream &out, const BernoulliStartState<DOMAIN> &startState) {
        out << "{ ";
        for(const auto &constraint: startState.constraints) {
            out << "P(X" << constraint.coefficients.indices[0] << ")=" << constraint.constant << " ";
        }
        for(const auto &factor: startState.logFactors) {
            out << "P(X" << factor.dependencies[0] << ")=" << 1.0 - factor(DOMAIN::zero).first << " ";
        }
        out << "}";
        return out;
    }

//private:
//    friend class boost::serialization::access;
//
//    template <typename Archive>
//    void load(Archive &ar, const unsigned int version) {
//        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//            double p;
//            ar >> p;
//            addBernoulliFactor(agentId, p);
//        }
//    }
//
//    template <typename Archive>
//    void save(Archive &ar, const unsigned int version) const {
//        for(int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
//            ar << probability(agentId);
//        }
//    }
//
//    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

#endif //ABMCMC_BERNOULLISTARTSTATE_H
