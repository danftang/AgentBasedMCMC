//
// Created by daniel on 18/08/2021.
//

#ifndef GLPKTEST_TRAJECTORYLIKELIHOODPMF_H
#define GLPKTEST_TRAJECTORYLIKELIHOODPMF_H

#include "ConvexPMF.h"
#include "AgentStateObservation.h"

//template<typename AGENT>
//class TrajectoryLikelihoodPMF {
//public:
//
//    TrajectoryLikelihoodPMF(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent)
//    : ConvexPMF(generateObservationLikelihood(realTrajectory, pMakeObservation, pObserveIfPresent)) { }
//
//    TrajectoryLikelihoodPMF(int nTimesteps, const AgentStateObservation<AGENT> &observation)
//    : ConvexPMF(likelihood(nTimesteps, observation)) { }
//
//    TrajectoryLikelihoodPMF(int nTimesteps)
//    : ConvexPMF([](const std::vector<double> &X) { return 1.0; }, Trajectory<AGENT>::dimension(nTimesteps)) { }
//
//
//    static ConvexPMF generateObservationLikelihood(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
//        int nTimesteps = realTrajectory.nTimesteps();
//        ConvexPMFProduct observations(realTrajectory.size());
//        for (int t=0; t<nTimesteps;++t) {
//            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
//                if (Random::nextDouble() < pMakeObservation) {
//                    AGENT agent(agentId);
//                    int nObserved = Random::nextBinomial(realTrajectory(t,agent), pObserveIfPresent);
//                    AgentStateObservation<AGENT> observation(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
//                    observations *= ConvexPMF<AGENT>(nTimesteps, observation);
//                }
//            }
//        }
//        //        checkTrajectorySatisfiesObervations(exactEndState, observations)
//        return observations;
//    }
//
//    static ConvexPMF likelihood(int nTimesteps, const AgentStateObservation<AGENT> &observation) {
//        ConvexPMF::PMF logP;
//        if(observation.state.time == nTimesteps) {
//            logP = [observation](const std::vector<double> &X) {
//                return observation.logP(observation.state.backwardOccupationNumber(X));
//            };
//        } else {
//            logP = [observation](const std::vector<double> &X) {
//                return observation.logP(observation.state.forwardOccupationNumber(X));
//            };
//        }
//        return ConvexPMF(std::move(logP), Trajectory<AGENT>::dimension(nTimesteps), observation.support());
//    }
//
//
//};


#endif //GLPKTEST_TRAJECTORYLIKELIHOODPMF_H
