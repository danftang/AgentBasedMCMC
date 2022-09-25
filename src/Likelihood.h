// Represents the likelihood function for a set of noisy and noiseless observations of
// agent state. Thw noisy observations are stored as factors in the distribution, and
// in noisyObservations (for serialisation), while the noiseless observations are
// stored in the constraints of the distribution.
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_LIKELIHOOD_H
#define ABMCMC_LIKELIHOOD_H

#include <boost/math/distributions/binomial.hpp>
#include "ConstrainedFactorisedDistribution.h"

template<class DOMAIN, class AGENT = typename DOMAIN::agent_type>
class Likelihood: public ConstrainedFactorisedDistribution<DOMAIN> {
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
    Likelihood(const DOMAIN &realTrajectory, double pMakeObservation, double pObserveIfPresent):
            Likelihood(generateObservations(realTrajectory, pMakeObservation, pObserveIfPresent), pObserveIfPresent)
    { }


    // single observation
    Likelihood(const State<AGENT> &state, typename DOMAIN::value_type nObserved, double pObserveIfPresent) {
        addObservation(state, nObserved, pObserveIfPresent);
    }


    // binomial (n k) p^k (1-p)^(n-k)
    // logBinomial klog(p) + (n-k)log(1-p) + log(n choose k)
    void addObservation(const State<AGENT> &state, typename DOMAIN::value_type nObserved, double pObserveIfPresent) {
        double logPnObserved = nObserved*log(pObserveIfPresent);
        if(pObserveIfPresent == 1.0) {
            this->constraints.push_back(X(DOMAIN::indexOf(state)) == nObserved); // noiseless
        } else {
            this->addFactor(
                    [state, nObserved, pObserveIfPresent, logPnObserved](const DOMAIN &trajectory) {
                        typename DOMAIN::value_type realOccupation = trajectory[state];
                        if (realOccupation < nObserved)
                            return std::pair(logPnObserved + ABM::kappa * (realOccupation - nObserved), false);
                        return std::pair(log(boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent),
                                                              nObserved)), true);
                    },
                    { DOMAIN::indexOf(state) }
            );
        }
    }

    static std::vector<std::pair<State<AGENT>,int>> generateObservations(const DOMAIN &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        std::vector<std::pair<State<AGENT>,int>> observations;
        for (int t = 0; t < realTrajectory.nTimesteps; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize; ++agentId) {
                if (Random::nextBool(pMakeObservation)) {
                    State<AGENT> state(t, agentId);
                    int nObserved = Random::nextBinomial(realTrajectory[state], pObserveIfPresent);
                    observations.emplace_back(state, nObserved);
                }
            }
        }
        return observations;
    }
};


#endif //ABMCMC_LIKELIHOOD_H
