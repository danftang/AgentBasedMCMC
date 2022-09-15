//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_BERNOULLISTARTSTATE_H
#define ABMCMC_BERNOULLISTARTSTATE_H

#include "ConstrainedFactorisedDistribution.h"

template<class AGENT>
class BernoulliStartState: public ConstrainedFactorisedDistribution<ModelState<AGENT>> {
public:

    explicit BernoulliStartState(std::function<double(AGENT)> agentToProb) {
        this->logFactors.reserve(AGENT::domainSize());
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            addBernoulliFactor(agentId, agentToProb(AGENT(agentId)));
        }
    }

    explicit BernoulliStartState(std::initializer_list<double> probs): BernoulliStartState([&probs](AGENT agent) {
        return data(probs)[agent];
    }) {}


    std::function<const ModelState<AGENT> &()> sampler() const {
        return [this, sample = ModelState<AGENT>()]() mutable -> const ModelState<AGENT> & {
            for(const auto &constraint: this->constraints)
                sample[constraint.coefficients.indices[0]] = constraint.constant;

            for(const auto &factor: this->logFactors) {
                double p = 1.0 - exp(factor(ModelState<AGENT>::zero).first);
                sample[factor.dependencies[0]] = Random::nextBool(p) ? 1:0;
//                std::cout << "p = " << p << " " << factor.dependencies[0] << " -> " << val << std::endl;
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
        if(p != 1.0 && p != 0.0) {
            this->addFactor(
                    SparseFunction<std::pair<double, bool>, const ModelState<AGENT> &>(
                            [p, agentId](const ModelState<AGENT> &modelState) {
                                return BernoulliStartState<AGENT>
                                ::widenedLogBernoulli(p, modelState[agentId]); // log of Bernoulli
                            },
                            {agentId}
                    )
            );
        } else {
            this->addConstraint(1*X(agentId) == (int)p);
        }
    }

//    double probability(int agentId) const {
//        return 1.0 - exp(this->logFactors[agentId](ModelState<AGENT>::zero).first);
//    }

    // convert to a distribution over trajectories
    operator ConstrainedFactorisedDistribution<Trajectory<AGENT>>() {
        return this->toTrajectoryDistribution(0);
    }


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
    static std::pair<double,bool> widenedLogBernoulli(double p, int occupation) {
        if(occupation == 0) return std::pair(log(1-p), true);
        if(occupation < 0) return std::pair(ABM::kappa*occupation, false);
        if(occupation > 1) return std::pair(ABM::kappa*(1-occupation), false);
        return std::pair(log(p), true);
    }

    friend std::ostream &operator <<(std::ostream &out, const BernoulliStartState<AGENT> &startState) {
        out << "{ ";
        for(const EqualityConstraint<ABM::occupation_type> &constraint: startState.constraints) {
            out << "P(X" << constraint.coefficients.indices[0] << ")=" << constraint.constant << " ";
        }
        for(const auto &factor: startState.logFactors) {
            out << "P(X" << factor.dependencies[0] << ")=" << 1.0 - exp(factor(ModelState<AGENT>::zero).first) << " ";
        }
        out << "}";
        return out;
    }

//private:
//    friend class boost::serialization::access;
//
//    template <typename Archive>
//    void load(Archive &ar, const unsigned int version) {
//        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//            double p;
//            ar >> p;
//            addBernoulliFactor(agentId, p);
//        }
//    }
//
//    template <typename Archive>
//    void save(Archive &ar, const unsigned int version) const {
//        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//            ar << probability(agentId);
//        }
//    }
//
//    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

#endif //ABMCMC_BERNOULLISTARTSTATE_H
