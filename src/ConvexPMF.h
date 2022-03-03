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
#include "ABMConstraints.h"

template<typename T> class ConvexPMFProduct;

// Represents a probability mass function defined on the vertices of a convex polyhedron.
template<typename DOMAIN>
class ConvexPMF: public ConvexPMFBase<DOMAIN> {
public:
    ConvexPMF(
            std::function<double(const DOMAIN &)> extendedLogP,
            int nDimensions,
            ConvexPolyhedron constraints = ConvexPolyhedron()
    )
    : ConvexPMFBase<DOMAIN>(extendedLogP, nDimensions, constraints) {};

    ConvexPMF(const ConvexPMFBase<DOMAIN> &base): ConvexPMFBase<DOMAIN>(base) { }
    ConvexPMF(ConvexPMFBase<DOMAIN> &&base): ConvexPMFBase<DOMAIN>(std::move(base)) { }

};


template<typename AGENT>
class ConvexPMF<Trajectory<AGENT>>: public ConvexPMFBase<Trajectory<AGENT>> {
public:
    ConvexPMF(
        std::function<double(const Trajectory<AGENT> &)> extendedLogP,
        int nDimensions,
        ConvexPolyhedron constraints = ConvexPolyhedron()
    )
    : ConvexPMFBase<Trajectory<AGENT>>(extendedLogP, nDimensions, constraints) {
        assert((nDimensions-1)%(AGENT::domainSize()*AGENT::actDomainSize()) == 0);
    };

    ConvexPMF(const ConvexPMFBase<Trajectory<AGENT>> &base): ConvexPMFBase<Trajectory<AGENT>>(base) { }
    ConvexPMF(ConvexPMFBase<Trajectory<AGENT>> &&base): ConvexPMFBase<Trajectory<AGENT>>(std::move(base)) { }

    ConvexPMF(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent)
    : ConvexPMF(generateObservationLikelihood(realTrajectory, pMakeObservation, pObserveIfPresent)) { }


    ConvexPMF(int nTimesteps, const AgentStateObservation<AGENT> &observation)
    : ConvexPMF(likelihood(nTimesteps, observation)) { }

    ConvexPMF(int nTimesteps, std::vector<AgentStateObservation<AGENT>> observations)
            : ConvexPMF(likelihood(nTimesteps, std::move(observations))) { }


    ConvexPMF(int nTimesteps, ConvexPMF<ModelState<AGENT>> priorStartState)
    : ConvexPMF([startState = std::move(priorStartState.extendedLogProb)](const Trajectory<AGENT> &T) {
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
//        int nTimesteps = realTrajectory.nTimesteps();
//        ConvexPMFProduct<Trajectory<AGENT>> observations(realTrajectory.size());
//        for (int t=0; t<nTimesteps;++t) {
//            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
//                if (Random::nextDouble() < pMakeObservation) {
//                    AGENT agent(agentId);
//                    int nObserved = Random::nextBinomial(realTrajectory(t,agent), pObserveIfPresent);
//                    AgentStateObservation<AGENT> observation(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
////                    debug(std::cout << "Adding constraint " << observation.support() << std::endl);
//                    observations *= likelihood(nTimesteps, observation);
//                }
//            }
//        }
//        assert(observations.convexSupport.isValidSolution(realTrajectory));
//        //        checkTrajectorySatisfiesObervations(exactEndState, observations)
//        return observations;
        return ConvexPMF<Trajectory<AGENT>>(
                realTrajectory.nTimesteps(),
                AgentStateObservation<AGENT>::generateObservations(
                        realTrajectory,
                        pMakeObservation,
                        pObserveIfPresent
                )
        );
    }


    static ConvexPMF<Trajectory<AGENT>> likelihood(int nTimesteps, const AgentStateObservation<AGENT> &observation) {
        std::function<double(const Trajectory<AGENT> &)> extendedLogP;
        if(observation.state.time == nTimesteps) {
            extendedLogP = [observation](const Trajectory<AGENT> &X) {
                return observation.extendedLogP(observation.state.backwardOccupationNumber(X));
            };
        } else {
            extendedLogP = [observation](const Trajectory<AGENT> &X) {
                return observation.extendedLogP(observation.state.forwardOccupationNumber(X));
            };
        }
        return ConvexPMF(std::move(extendedLogP), Trajectory<AGENT>::dimension(nTimesteps), observation.support());
    }

    static ConvexPMF<Trajectory<AGENT>> likelihood(int nTimesteps, std::vector<AgentStateObservation<AGENT>> observations) {
        ConvexPolyhedron support;
        for(const AgentStateObservation<AGENT> &observation: observations) {
            support += observation.support();
        }
        std::function<double(const Trajectory<AGENT> &)> extendedLogP = [observations = std::move(observations)](const Trajectory<AGENT> &X) {
            double extLogP = 0.0;
            for(const AgentStateObservation<AGENT> &observation: observations) {
                extLogP += observation.extendedLogP(observation.state.occupationNumber(X));
            }
            return extLogP;
        };
        return ConvexPMF(std::move(extendedLogP), Trajectory<AGENT>::dimension(nTimesteps), std::move(support));
    }


};


template<typename AGENT>
class ConvexPMF<ModelState<AGENT>>: public ConvexPMFBase<ModelState<AGENT>> {
public:
    ConvexPMF(std::function<double(const ModelState<AGENT> &)> extendedLogP, ConvexPolyhedron constraints = ConvexPolyhedron())
    : ConvexPMFBase<ModelState<AGENT>>(extendedLogP, AGENT::domainSize(), constraints) {};


    ConvexPMF(std::function<double(const ModelState<AGENT> &)> extendedLogP, int nDimensions, ConvexPolyhedron constraints = ConvexPolyhedron())
    : ConvexPMF(extendedLogP, constraints) {
        assert(nDimensions == AGENT::domainSize());
    };

    ConvexPMF(const ConvexPMFBase<ModelState<AGENT>> &base): ConvexPMFBase<ModelState<AGENT>>(base) { }
    ConvexPMF(ConvexPMFBase<ModelState<AGENT>> &&base): ConvexPMFBase<ModelState<AGENT>>(std::move(base)) { }

//    ConvexPMF(const std::function<double(AGENT, int)> &marginalLogProb)
//    : ConvexPMF(
//        [marginalLogProb](const ModelState<AGENT> &M) {
//            double extendedLogP = 0.0;
//            for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
//                extendedLogP += extendedMarginalLogProb(AGENT(agentId),std::round(M[agentId]));
//            }
//            return extendedLogP;
//        },
//        constraints(marginalLogProb)
//    ) {}
//
//
//    static ConvexPolyhedron constraints(const std::function<double(AGENT, int)> &marginalLogProb) {
//        ConvexPolyhedron constraints;
//        for(int agentId=0; agentId<AGENT::domainSize(); ++agentId) {
//            int lowerBound = 0;
//            while(marginalLogProb(AGENT(agentId), lowerBound) <= -DBL_MAX) ++lowerBound;
//            int upperBound = lowerBound;
//            double cumulativeP = exp(marginalLogProb(AGENT(agentId), lowerBound));
//            while(cumulativeP < 1.0-1e-6) {
//                cumulativeP += exp(marginalLogProb(AGENT(agentId), ++upperBound));
//            }
//            constraints.push_back(lowerBound <= 1.0*glp::X(agentId) <= upperBound);
//        }
//        return constraints;
//    }

};



#endif //GLPKTEST_CONVEXPMF_H
