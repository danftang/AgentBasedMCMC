// Represents the likelihood function for a set of noisy and noiseless observations of
// agent state. Thw noisy observations are stored as factors in the distribution, and
// in noisyObservations (for serialisation), while the noiseless observations are
// stored in the constraints of the distribution.
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_LIKELIHOOD_H
#define ABMCMC_LIKELIHOOD_H

#include "ABM.h"
#include "NoisyAgentStateObservation.h"
#include "NoiselessAgentStateObservation.h"
#include "Trajectory.h"
#include "ConstrainedFactorisedDistribution.h"

template<class AGENT>
class Likelihood: public ConstrainedFactorisedDistribution<Trajectory<AGENT>,ABM::coefficient_type> {
public:
    Likelihood()=default;

    Likelihood(const std::vector<std::pair<State<AGENT>,int>> &observations, double pObserveIfPresent) {
        assert(pObserveIfPresent != 0.0); // pointless as no information
        if(pObserveIfPresent != 1.0) {
            this->factors.reserve(observations.size());
        } else {
            this->constraints.reserve(observations.size());
        }
        for(const std::pair<State<AGENT>,int> &observation: observations) {
            addObservation(observation.first, observation.second, pObserveIfPresent);
        }
    }

    // generate random observations
    Likelihood(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent):
            Likelihood(generateObservations(realTrajectory, pMakeObservation, pObserveIfPresent), pObserveIfPresent)
    { }


    // single observation
    Likelihood(const State<AGENT> &state, ABM::occupation_type nObserved, double pObserveIfPresent) {
        addObservation(state, nObserved, pObserveIfPresent);
    }


    // binomial (n k) p^k (1-p)^(n-k)
    // logBinomial klog(p) + (n-k)log(1-p) + log(n choose k)
    void addObservation(const State<AGENT> &state, ABM::occupation_type nObserved, double pObserveIfPresent) {
        double logPnObserved = nObserved*log(pObserveIfPresent);
        if(pObserveIfPresent == 1.0) {
            this->addConstraint(1*state == nObserved); // noiseless
        } else {
            this->addFactor(
                    SparseFunction<std::pair<double,bool>, const Trajectory<AGENT> &>(
                            [state, nObserved, pObserveIfPresent, logPnObserved](const Trajectory<AGENT> &trajectory) {
                                ABM::occupation_type realOccupation = trajectory[state];
                                if(realOccupation < nObserved) return std::pair(logPnObserved + ABM::kappa*(realOccupation - nObserved),false);
                                return std::pair(log(boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), nObserved)),true);
                            },
                            state.forwardOccupationDependencies()
                            )
                    );
        }
    }

    static std::vector<std::pair<State<AGENT>,int>> generateObservations(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        std::vector<std::pair<State<AGENT>,int>> observations;
        for (int t = 0; t < realTrajectory.nTimesteps(); ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextBool(pMakeObservation)) {
                    State<AGENT> state(t, agentId);
                    int nObserved = Random::nextBinomial(realTrajectory[state], pObserveIfPresent);
                    observations.emplace_back(state, nObserved);
                }
            }
        }
        return observations;
    }

//    friend std::ostream &operator <<(std::ostream &out, const Likelihood<AGENT> &likelihood) {
//        for(const EqualityConstraint<ABM::occupation_type> &noiselessObservation : likelihood.constraints) {
//            out << noiselessObservation << std::endl;
//        }
//        for(const NoisyAgentStateObservation<AGENT> &observation : likelihood.noisyObservations) {
//            out << observation << std::endl;
//        }
//        return out;
//    }

//private:
//    friend class boost::serialization::access;
//
//    template <typename Archive>
//    void load(Archive &ar, const unsigned int version) {
//        ar >> noisyObservations >> this->constraints;
//        for(const auto &observation: noisyObservations) this->addFactor(observation.toSparseWidenedFunction());
//    }
//
//    template <typename Archive>
//    void save(Archive &ar, const unsigned int version) const {
//        ar << noisyObservations << this->constraints;
//    }
//
//    BOOST_SERIALIZATION_SPLIT_MEMBER();
};


#endif //ABMCMC_LIKELIHOOD_H
