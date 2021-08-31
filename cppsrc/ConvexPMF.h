//
// Created by daniel on 26/07/2021.
//

#ifndef GLPKTEST_CONVEXPMF_H
#define GLPKTEST_CONVEXPMF_H

#include "glpkpp.h"
#include "ConvexPolyhedron.h"
#include "ConvexPMFBase.h"
#include "Trajectory.h"
#include "ConvexPMFProduct.h"
#include "AgentStateObservation.h"

template<typename T> class ConvexPMFProduct;

// Represents a probability mass function defined on the vertices of a convex polyhedron.
template<typename DOMAIN>
class ConvexPMF: public ConvexPMFBase<DOMAIN> {
public:
    ConvexPMF(
            std::function<double(const DOMAIN &)> logP,
            int nDimensions,
            ConvexPolyhedron constraints = ConvexPolyhedron()
    )
    : ConvexPMFBase<DOMAIN>(logP, nDimensions, constraints) {};

    ConvexPMF(const ConvexPMFBase<DOMAIN> &base): ConvexPMFBase<DOMAIN>(base) { }
    ConvexPMF(ConvexPMFBase<DOMAIN> &&base): ConvexPMFBase<DOMAIN>(std::move(base)) { }

};


template<typename AGENT>
class ConvexPMF<Trajectory<AGENT>>: public ConvexPMFBase<Trajectory<AGENT>> {
public:
    ConvexPMF(
        std::function<double(const Trajectory<AGENT> &)> logP,
        int nDimensions,
        ConvexPolyhedron constraints = ConvexPolyhedron()
    )
    : ConvexPMFBase<Trajectory<AGENT>>(logP, nDimensions, constraints) {};

    ConvexPMF(const ConvexPMFBase<Trajectory<AGENT>> &base): ConvexPMFBase<Trajectory<AGENT>>(base) { }
    ConvexPMF(ConvexPMFBase<Trajectory<AGENT>> &&base): ConvexPMFBase<Trajectory<AGENT>>(std::move(base)) { }

    ConvexPMF(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent)
    : ConvexPMF(generateObservationLikelihood(realTrajectory, pMakeObservation, pObserveIfPresent)) { }


    ConvexPMF(int nTimesteps, const AgentStateObservation<AGENT> &observation)
    : ConvexPMF(likelihood(nTimesteps, observation)) { }


    ConvexPMF(int nTimesteps, ConvexPMF<ModelState<AGENT>> priorStartState)
    : ConvexPMF([startState = std::move(priorStartState.logProb)](const Trajectory<AGENT> &T) {
        return T.logProb() + startState(T(0));
        },
                Trajectory<AGENT>::dimension(nTimesteps),
                ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps) +
                ABMConstraints<AGENT>::startStateConstraintsToTrajectoryConstraints(priorStartState.convexSupport)) {}



    static ConvexPMF<Trajectory<AGENT>> uniformTrajectoryDistribution(int nTimesteps) {
        return ConvexPMF([](const Trajectory<AGENT> &T) { return 1.0; }, Trajectory<AGENT>::dimension(nTimesteps));
    }


    static ConvexPMF<Trajectory<AGENT>>
    generateObservationLikelihood(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        int nTimesteps = realTrajectory.nTimesteps();
        ConvexPMFProduct<Trajectory<AGENT>> observations(realTrajectory.size());
        for (int t=0; t<nTimesteps;++t) {
            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    AGENT agent(agentId);
                    int nObserved = Random::nextBinomial(realTrajectory(t,agent), pObserveIfPresent);
                    AgentStateObservation<AGENT> observation(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
                    observations *= likelihood(nTimesteps, observation);
                }
            }
        }
        //        checkTrajectorySatisfiesObervations(exactEndState, observations)
        return observations;
    }


    static ConvexPMF<Trajectory<AGENT>> likelihood(int nTimesteps, const AgentStateObservation<AGENT> &observation) {
        std::function<double(const Trajectory<AGENT> &)> logP;
        if(observation.state.time == nTimesteps) {
            logP = [observation](const Trajectory<AGENT> &X) {
                return observation.logP(observation.state.backwardOccupationNumber(X));
            };
        } else {
            logP = [observation](const Trajectory<AGENT> &X) {
                return observation.logP(observation.state.forwardOccupationNumber(X));
            };
        }
        return ConvexPMF(std::move(logP), Trajectory<AGENT>::dimension(nTimesteps), observation.support());
    }


};


template<typename AGENT>
class ConvexPMF<ModelState<AGENT>>: public ConvexPMFBase<ModelState<AGENT>> {
public:
    ConvexPMF(std::function<double(const ModelState<AGENT> &)> logP, ConvexPolyhedron constraints = ConvexPolyhedron())
    : ConvexPMFBase<ModelState<AGENT>>(logP, AGENT::domainSize(), constraints) {};


    ConvexPMF(std::function<double(const ModelState<AGENT> &)> logP, int nDimensions, ConvexPolyhedron constraints = ConvexPolyhedron())
    : ConvexPMF(logP, constraints) {
        assert(nDimensions == AGENT::domainSize());
    };


    ConvexPMF(const ConvexPMFBase<ModelState<AGENT>> &base): ConvexPMFBase<ModelState<AGENT>>(base) { }
    ConvexPMF(ConvexPMFBase<ModelState<AGENT>> &&base): ConvexPMFBase<ModelState<AGENT>>(std::move(base)) { }

    ConvexPMF(const std::function<double(AGENT, int)> &marginalLogProb)
    : ConvexPMF(
        [marginalLogProb](const ModelState<AGENT> &M) {
            double logP = 0.0;
            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                logP += marginalLogProb(AGENT(agentId),std::round(M[agentId]));
            }
            return logP;
        },
        constraints(marginalLogProb)
    ) {}


    static ConvexPolyhedron constraints(const std::function<double(AGENT, int)> &marginalLogProb) {
        ConvexPolyhedron constraints;
        for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
            int lowerBound = 0;
            while(marginalLogProb(AGENT(agentId), lowerBound) <= -DBL_MAX) ++lowerBound;
            int upperBound = lowerBound;
            double cumulativeP = exp(marginalLogProb(AGENT(agentId), lowerBound));
            while(cumulativeP < 1.0-1e-6) {
                cumulativeP += exp(marginalLogProb(AGENT(agentId), ++upperBound));
            }
            constraints.push_back(lowerBound <= 1.0*glp::X(agentId) <= upperBound);
        }
        return constraints;
    }

};



#endif //GLPKTEST_CONVEXPMF_H
