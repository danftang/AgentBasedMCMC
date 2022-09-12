//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_BERNOULLISTARTSTATE_H
#define ABMCMC_BERNOULLISTARTSTATE_H

#include "State.h"
#include "StartStateDistribution.h"

template<class AGENT>
class BernoulliStartState: public FactorisedDistribution<ModelState<AGENT>> {
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


    ModelState<AGENT> nextSample() const {
        ModelState<AGENT> sample;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            sample[agentId] = Random::nextDouble(probability(agentId))<probability(agentId)?1:0;
        }
        return sample;
    }


    // A Poisson distribution that decays exponentially below zero
    static std::pair<double,bool> widenedLogBernoulli(double p, int occupation) {
        if(occupation == 0) return std::pair(log(1-p), true);
        if(occupation < 0) return std::pair(ABM::kappa*occupation, false);
        if(occupation > 1) return std::pair(ABM::kappa*(1-occupation), false);
        return std::pair(log(p), true);
    }


    using FactorisedDistribution<ModelState<AGENT>>::addFactor;

    void addBernoulliFactor(int agentId, double p) {
        addFactor(
                SparseFunction<std::pair<double, bool>, const ModelState<AGENT> &>(
                        [p, agentId](const ModelState<AGENT> &modelState) {
                            return BernoulliStartState<AGENT>
                                   ::widenedLogBernoulli(p, modelState[agentId]); // log of Bernoulli
                        },
                        {agentId}
                )
        );
    }

    double probability(int agentId) const {
        return 1.0 - exp(this->logFactors[agentId](ModelState<AGENT>::zero).first);
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
        double p = probability(agentId);
        return SparseFunction<std::pair<double,bool>,const Trajectory<AGENT> &>(
                [p,agentId](const Trajectory<AGENT> &trajectory) {
                    return BernoulliStartState<AGENT>::widenedLogBernoulli(p, trajectory[State<AGENT>(0,agentId)]);
                },
                State<AGENT>(0,agentId).forwardOccupationDependencies()
        );
    }


    SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &> getModelStateFactor(int agentId) {
        double p = probability(agentId);
        return SparseFunction<std::pair<double,bool>,const ModelState<AGENT> &>(
                [p,agentId](const ModelState<AGENT> &modelState) {
                    return BernoulliStartState<AGENT>::widenedLogBernoulli(p, modelState[agentId]); // log of Poisson
                },
                {agentId}
        );
    }

    friend std::ostream &operator <<(std::ostream &out, const BernoulliStartState<AGENT> &startState) {
        out << "{ ";
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            out << startState.probability(agentId) << " ";
        }
        out << "}";
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            double p;
            ar >> p;
            addBernoulliFactor(agentId, p);
        }
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            ar << probability(agentId);
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};

#endif //ABMCMC_BERNOULLISTARTSTATE_H
