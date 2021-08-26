//
// Created by daniel on 27/07/2021.
//

#ifndef GLPKTEST_TRAJECTORYPRIORDISTRIBUTION_H
#define GLPKTEST_TRAJECTORYPRIORDISTRIBUTION_H

#include <math.h>
#include "ConvexPMF.h"
#include "Trajectory.h"
#include "ABMConstraints.h"
#include "AgentStateObservation.h"
#include "Distribution.h"
#include "SampleStatistics.h"

// Represents the (measurable) space of trajectories of an ABM with a given number of timesteps,
// from which we can generate various useful PMFs
//template<typename AGENT>
//class TrajectoryPriorDistribution: public Distribution {
//public:
//    int             nTimesteps;
//    ConvexPMF       startStatePMF;
//    std::function<std::vector<double>()> startStateSampler;
//
//    TrajectoryPriorDistribution(int nTimesteps, const Distribution &startStateDistribution)
//    : nTimesteps(nTimesteps),
//    startStatePMF(startStateDistribution.PMF()),
//    startStateSampler(startStateDistribution.sampler())
//    {
//        assert(nTimesteps != 0); // don't allow null space
//    }
//
//
//    // Takes a PMF over the start state of an ABM and returns the PMF over trajectories of this window length
//    // such that the probability of a solution is the trajectoryPrior probability of the solution times the probability
//    // of the start state.
//    ConvexPMF PMF() const {
//        return ConvexPMF([startState = startStatePMF.logProb](const std::vector<double> &X) {
//            const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
//            return T.logProb() + startState(T(0));
//            },
//                         nDimensions(),
//                         ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps) +
//                         ABMConstraints<AGENT>::startStateConstraintsToTrajectoryConstraints(startStatePMF.convexSupport));
//    }
//
//
//    std::function<std::vector<double>()> sampler() const {
//        return [startSampler = startStateSampler,nTimesteps = this->nTimesteps]() {
//            return Trajectory<AGENT>(nTimesteps, startSampler);
//        };
//    }
//
//
////    ConvexPMF actFermmionicDistribution() {
////        return ConvexPMF(
////                Trajectory<AGENT>::logProb,
////                nDimensions(),
////                ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps));
////    }
//
//
//    // turns a number of 1-dimensional AGENT state distributions into a trajectory distribution
//    static ConvexPMF startStateDistribution(std::function<ConvexPMF(AGENT)> agentStateDistribution) {
//
//    }
//
//
//
////    ConvexPMF likelihood(const AgentStateObservation<AGENT> &observation) const {
////        ConvexPMF::PMF logP;
////        if(observation.state.time == nTimesteps) {
////            logP = [observation](const std::vector<double> &X) {
////                return observation.logP(observation.state.backwardOccupationNumber(X));
////            };
////        } else {
////            logP = [observation](const std::vector<double> &X) {
////                return observation.logP(observation.state.forwardOccupationNumber(X));
////            };
////        }
////        return ConvexPMF(std::move(logP), nDimensions(), observation.support());
////    }
//
//
//    int nDimensions() const {
//        return nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()+1;
//    }
//
////    static ConvexPMF generateObservationLikelihood(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
////        int nTimesteps = realTrajectory.nTimesteps();
////        ConvexPMFProduct observations(realTrajectory.size());
////        TrajectoryPriorDistribution<AGENT> window(realTrajectory.nTimesteps());
////        for (int t=0; t<nTimesteps;++t) {
////            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
////                if (Random::nextDouble() < pMakeObservation) {
////                    AGENT agent(agentId);
////                    int nObserved = Random::nextBinomial(realTrajectory(t,agent), pObserveIfPresent);
////                    AgentStateObservation<AGENT> observation(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
////                    observations *= window.likelihood(observation);
////                }
////            }
////        }
////        //        checkTrajectorySatisfiesObervations(solution, observations)
////        return observations;
////    }
//
//
//    // Fit a Binomial to each state such that n = the max observed occupation
//    // and the mean = np = the observed mean.
////    std::vector<double> priorEndStateStats(int nSamples) {
////        SampleStatistics endStateStats;
////        for(int s=0; s<nSamples; ++s) {
////            Trajectory<AGENT> trajectory(nTimesteps(), startStatePrior.nextSample());
////            int t = 0;
////            for(int w=0; w<windows.size(); ++w) {
////                t += windows[w].realTrajectory.nTimesteps();
////                endOfWindowStates[w] += trajectory(t);
////            }
////        }
////        return endOfWindowStates;
////    }
//
//    // Calculates the information gained by assimilation over nTimesteps
//    //
//
////    double calculateInformationGain() {
////        std::vector<PoissonState<AGENT>> referencePriors = priorEndState(10000);
////        double gain = informationGain(
////                    windows[w].realTrajectory.endState(),
////                    referencePriors[w],
////                    windows[w].analysis);
////            Gnuplot gp;
////            Experiments::plotHeatMap(gp, referencePriors[w], windows[w].realTrajectory.endState());
////            Gnuplot gp2;
////            Experiments::plotHeatMap(gp2, windows[w].analysis, windows[w].realTrajectory.endState());
////        return gain;
////    }
//
//    // Information gain, in bits, about a real vector given a trajectoryPrior PMF and a posterior PMF
//    static double informationGain(const std::vector<double> &realState, const ConvexPMF &prior, const ConvexPMF &posterior) {
//        return (posterior.logProb(realState) - prior.logProb(realState))/log(2);
//    }
//
//
//
//};
//
//
#endif //GLPKTEST_TRAJECTORYPRIORDISTRIBUTION_H

