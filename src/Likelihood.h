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
    std::vector<NoisyAgentStateObservation<AGENT>> noisyObservations;

    Likelihood()=default;

    // generate random observations
    Likelihood(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        int nTimesteps = realTrajectory.nTimesteps();
        int approxSize = nTimesteps * AGENT::domainSize() * pMakeObservation;
        if(pMakeObservation != 1.0) {
            this->logFactors.reserve(approxSize);
            noisyObservations.reserve(approxSize);
        } else {
            this->constraints.reserve(approxSize);
        }
        for (int t = 0; t < nTimesteps; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    State<AGENT> state(t, agentId);
                    int nObserved = Random::nextBinomial(realTrajectory[state], pObserveIfPresent);
                    addObservation(state, nObserved, pObserveIfPresent);
                }
            }
        }
    }


    void addObservation(const State<AGENT> &state, ABM::occupation_type nObserved, double pObserveIfPresent) {
        if(pObserveIfPresent == 1.0) {
            this->addConstraint(1*state == nObserved); // noiseless
        } else {
            noisyObservations.push_back(NoisyAgentStateObservation<AGENT>(state, nObserved, pObserveIfPresent));
            this->addFactor(noisyObservations.back().toSparseWidenedFunction());
        }
    }


    friend std::ostream &operator <<(std::ostream &out, const Likelihood<AGENT> &likelihood) {
        for(const EqualityConstraint<ABM::occupation_type> &noiselessObservation : likelihood.constraints) {
            out << noiselessObservation << std::endl;
        }
        for(const NoisyAgentStateObservation<AGENT> &observation : likelihood.noisyObservations) {
            out << observation << std::endl;
        }
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        ar >> noisyObservations >> this->constraints;
        for(const auto &observation: noisyObservations) this->addFactor(observation.toSparseWidenedFunction());
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar << noisyObservations << this->constraints;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();
};


#endif //ABMCMC_LIKELIHOOD_H
