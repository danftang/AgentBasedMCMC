//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_DATAASSIMILATION_H
#define GLPKTEST_DATAASSIMILATION_H

#include <vector>
#include "ABMProblem.h"
#include "Trajectory.h"
#include "PoissonState.h"
#include "AssimilationWindow.h"
#include "debug.h"

template <typename AGENT>
class DataAssimilation {
public:
    PoissonState<AGENT>                     startStatePrior;
    std::vector<AssimilationWindow<AGENT>>  windows;


//    DataAssimilation(const PoissonState<AGENT> &startStatePrior): startStatePrior(startStatePrior) {
//    }

    DataAssimilation(
            int nWindows,
            int nTimestepsPerWindow,
            const PoissonState<AGENT> &priorStartState,
            double pMakeObservation,
            double pObserveIfPresent,
            int nSamples): startStatePrior(priorStartState) {
        ModelState<AGENT> realState = startStatePrior.sample();
        const PoissonState<AGENT> *lastAnalysis = &startStatePrior;
        for(int t = 0; t < nWindows; ++t) {
            const AssimilationWindow<AGENT> &window =
                    windows.emplace_back(nTimestepsPerWindow, *lastAnalysis, realState, pMakeObservation, pObserveIfPresent, nSamples);
            lastAnalysis = &window.analysis;
            realState = window.realTrajectory(nTimestepsPerWindow);
            debug(        Gnuplot gp; Experiments::plotHeatMap(gp, poissonModelState, realState));
        }
    }



    // Calculates the information gained by assimilation over nTimesteps
    std::vector<double> calculateInformationGain() {
        std::vector<double> informationGains(windows.size());
        std::vector<PoissonState<AGENT>> referencePriors = priorWindows(AGENT::domainSize() * 50);
        for(int w=0; w<windows.size(); ++w) {
            informationGains[w] = informationGain(
                    windows[w].realTrajectory.endState(),
                    referencePriors[w],
                    windows[w].analysis);
        }
        return informationGains;
    }


//    static PoissonState<AGENT> assimilateWindow(
//            int nTimeseps,
//            std::vector<Observation<AGENT>> observations,
//            const PoissonState<AGENT> &startStatePrior,
////            const Trajectory<AGENT> &initialTrajectory,
//            int nSamples) {
//
//        ABMProblem<AGENT> abm(nTimeseps, observations, [&](const Trajectory<AGENT> &trajectory) {
//            return startStatePrior.logProb(trajectory(0));
//        });
//        SimplexMCMC mcmc(abm, abm.logProbFunc());
//        Trajectory<AGENT> initialTrajectory = generateInitialState(mcmc, startStatePrior);
//        mcmc.setLPState(initialTrajectory);
//        assert(mcmc.abmSanityChecks());
//
//        // test
//        std::vector<int> occupationHistogram(10,0.0);
//
//        PoissonState<AGENT> finalState;
//        for(int n=0; n<nSamples; ++n) {
//            mcmc.nextSample();
//            if(nSamples<1000 || n%1000 == 1) {
////                std::cout << "Sample " << n << std::endl;
//                assert(abm.isValidSolution(mcmc.X()));
////            std::cout << "Sample " << n << " : " << glp::SparseVec(mcmc.X()) << std::endl;
//            }
//            const Trajectory<AGENT> &trajectory = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
//            ModelState<AGENT> endStateSample = trajectory(abm.nTimesteps);
//            occupationHistogram[endStateSample[AGENT(152)]] += 1;
//            finalState += endStateSample;
//        }
//        debug(
//                std::cout << "Occupation histogram " << occupationHistogram << std::endl;
//                std::cout << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl;
//                std::cout << "Feasible sample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl;
//                std::cout << "Infeasible sample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
//        );
//        return finalState;
//    }


    // Calculates the posiion states at the end of each window, given no observations but
    // including the prior over the start state. This is needed in order to calculate the
    // information gain due to assimilation.
    std::vector<PoissonState<AGENT>> priorWindows(int nSamples) {
        int sumOfTimesteps = 0;
        for(const AssimilationWindow<AGENT> &window: windows) {
            sumOfTimesteps += window.realTrajectory.nTimesteps();
        }
        std::vector<PoissonState<AGENT>> endOfWindowStates(windows.size());
        for(int s=0; s<nSamples; ++s) {
            Trajectory<AGENT> trajectory(sumOfTimesteps, startStatePrior.sample());
            int t = 0;
            for(int w=0; w<windows.size(); ++w) {
                t += windows[w].realTrajectory.sumOfTimesteps();
                endOfWindowStates[w] += trajectory(t);
            }
        }
        return endOfWindowStates;
    }


//    static std::vector<PoissonState<AGENT>> posterior(
//            int nWindows,
//            int nTimestepsPerWindow,
//            int nSamplesPerWindow,
//            const PoissonState<AGENT> &startState,
//            const std::vector<std::vector<Observation<AGENT>>> &observations) {
//        ModelState<AGENT> realState = startState.sample();
//        std::vector<PoissonState<AGENT>> endOfWindowStates(nWindows);
//
//        const PoissonState<AGENT> *windowStart = &startState;
//        for (int window = 0; window < nWindows; ++window) {
//            endOfWindowStates[window] = DataAssimilation<AGENT>::assimilateWindow(
//                    nTimestepsPerWindow,
//                    observations,
//                    *windowStart,
//                    nSamplesPerWindow
//            );
//            windowStart = &endOfWindowStates.back();
//        }
//        return endOfWindowStates;
//    }


//    static std::vector<std::vector<Observation<AGENT>>> observe(
//            const Trajectory<AGENT> &realTrajectory,
//            int nWindows,
//            double pMakeObservation,
//            double pObserveIfPresent) {
//        std::vector<std::vector<Observation<AGENT>>> observations(nWindows);
//
//        const ModelState<AGENT> *windowStartState = &startState;
//        for (int window = 0; window < nWindows; ++window) {
//            auto obs = Observation<AGENT>::generateObservations(*windowStartState, nTimestepsPerWindow,
//                                                                                 pMakeObservation, pObserveIfPresent);
//            observations.push_back(obs);
//            endOfWindowStates.push_back(trajectory(nTimestepsPerWindow));
//            windowStartState = &endOfWindowStates.back();
//        }
//        return std::pair(observations, endOfWindowStates);
//    }


    static double  informationGain(const ModelState<AGENT> &realState, const PoissonState<AGENT> &prior, const PoissonState<AGENT> &posterior) {
        return (posterior.logProb(realState) - prior.logProb(realState))/log(2);
    }


};


#endif //GLPKTEST_DATAASSIMILATION_H
