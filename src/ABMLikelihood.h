// Represents the likelihood function for a set of noisy and noiseless observations of
// agent state. Thw noisy observations are stored as factors in the distribution, and
// in noisyObservations (for serialisation), while the noiseless observations are
// stored in the constraints of the distribution.
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_ABMLIKELIHOOD_H
#define ABMCMC_ABMLIKELIHOOD_H

#include <boost/math/distributions/binomial.hpp>
#include "ConstrainedFactorisedDistribution.h"

template<class TRAJECTORY, class AGENT = typename TRAJECTORY::agent_type>
class ABMLikelihood: public ConstrainedFactorisedDistribution<TRAJECTORY> {
public:
    std::vector<std::pair<State<AGENT>,int>> observations;
    double pObserveIfPresent;
    double kappa;
    std::optional<TRAJECTORY> realTrajectory; // just for the record

    ABMLikelihood()=default;

    ABMLikelihood(std::vector<std::pair<State<AGENT>,int>> observations, double pObserveIfPresent, double kappa):
        observations(std::move(observations)),
        pObserveIfPresent(pObserveIfPresent),
        kappa(kappa)
    {
        init();
    }

    // generate random observations
    ABMLikelihood(TRAJECTORY realTrajectory, double pMakeObservation, double pObserveIfPresent, double kappa):
        ABMLikelihood(generateObservations(realTrajectory, pMakeObservation, pObserveIfPresent), pObserveIfPresent, kappa)
    {
        realTrajectory = std::move(realTrajectory);
    }


    // single observation
    ABMLikelihood(const State<AGENT> &state, typename TRAJECTORY::value_type nObserved, double pObserveIfPresent, double kappa):
        ABMLikelihood({{state,nObserved}}, pObserveIfPresent, kappa)
    { }

    void init() {
        assert(pObserveIfPresent != 0.0); // pointless as no information
        if(pObserveIfPresent == 1.0) {
            this->constraints.reserve(observations.size());
        } else {
            this->factors.reserve(observations.size());
        }
        for(const std::pair<State<AGENT>,int> &observation: observations) {
            addObservation(observation.first, observation.second, pObserveIfPresent);
        }

    }

    // binomial (n k) p^k (1-p)^(n-k)
    // logBinomial klog(p) + (n-k)log(1-p) + log(n choose k)
    void addObservation(const State<AGENT> &state, typename TRAJECTORY::value_type nObserved, double pObserveIfPresent) {
        double logPnObserved = nObserved*log(pObserveIfPresent);
        if(pObserveIfPresent == 1.0) {
            this->constraints.push_back(X(TRAJECTORY::indexOf(state)) == nObserved); // noiseless
        } else {
            this->addFactor(
                    [state, nObserved, pObserveIfPresent, logPnObserved, kappa = this->kappa](const TRAJECTORY &trajectory) {
                        typename TRAJECTORY::value_type realOccupation = trajectory[state];
                        if(realOccupation < 0) std::pair(logPnObserved - kappa*nObserved, false);// TODO: TEST no decay in negative!!!!!
                        if (realOccupation < nObserved)
                            return std::pair(logPnObserved + kappa * (realOccupation - nObserved), false);
                        return std::pair(log(boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent),
                                                              nObserved)), true);
                    },
                    {TRAJECTORY::indexOf(state) }
            );
        }
    }

    static std::vector<std::pair<State<AGENT>,int>> generateObservations(const TRAJECTORY &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
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

    friend std::ostream &operator <<(std::ostream &out, const ABMLikelihood<TRAJECTORY> &likelihood) {
        out << "Observations: " << likelihood.observations << std::endl;
        out << "pObserveIfPresent: " << likelihood.pObserveIfPresent << std::endl;
        out << "kappa: " << likelihood.kappa << std::endl;
        out << "real trajectory: " << likelihood.realTrajectory << std::endl;
        out << static_cast<ConstrainedFactorisedDistribution<TRAJECTORY>>(likelihood);
        return out;
    }

private:

    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        bool realTrajectoryHasValue;
        ar >> observations >> pObserveIfPresent >> kappa >> realTrajectoryHasValue;
        if(realTrajectoryHasValue) {
            TRAJECTORY rTraj;
            ar >> rTraj;
            realTrajectory.template emplace(rTraj);
        }
        init();
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar << observations << pObserveIfPresent << kappa << realTrajectory.has_value();
        if(realTrajectory.has_value()) ar << realTrajectory.value();
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();


};


#endif //ABMCMC_ABMLIKELIHOOD_H
