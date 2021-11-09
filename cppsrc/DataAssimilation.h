//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_DATAASSIMILATION_H
#define GLPKTEST_DATAASSIMILATION_H

#include <vector>
#include "AssimilationWindow.h"
#include "debug.h"
#include "MCMCSolver.h"

//template <typename AGENT>
//class DataAssimilation {
//public:
//    std::vector<AssimilationWindow<AGENT>>  windows;
//    IntSampleStatistics analysis;
//    ConvexPMF analysisPMF;
//    std::function<std::vector<double>()> analysisSampler;
//    double pMakeObservation;
//    double pObserveIfPresent;
//
//    DataAssimilation(const Distribution &startStatePrior, double pMakeObservation, double pObserveIfPresent)
//    : analysis(AGENT::domainSize()),
//    analysisPMF(startStatePrior.PMF()),
//    analysisSampler(startStatePrior.sampler()),
//    pMakeObservation(pMakeObservation),
//    pObserveIfPresent(pObserveIfPresent)
//    {}
//
////
////    DataAssimilation(
////            PoissonState<AGENT> priorStartState,
////            std::function<std::vector<Observation<AGENT>>(const Trajectory<AGENT> &)> observationOperator):
////            startStatePrior(priorStartState),
////            observationOperator(observationOperator) { }
////
////
////    DataAssimilation(
////            int nWindows,
////            int nTimestepsPerWindow,
////            const PoissonState<AGENT> &priorStartState,
////            std::function<std::vector<Observation<AGENT>>(const Trajectory<AGENT> &)> observationOperator,
////            int nSamples,
////            int nBurnInSamples):
////            startStatePrior(priorStartState),
////            observationOperator(observationOperator)
////    {
////        for(int t = 0; t < nWindows; ++t) {
////            addWindow(nTimestepsPerWindow, nSamples, nBurnInSamples);
////            debug(        Gnuplot gp; Experiments::plotHeatMap(gp, windows[t].analysis, windows[t].realTrajectory.endState()));
////        }
////    }
////
////
//const AssimilationWindow<AGENT> &addWindow(int nTimesteps, int nBurnInSamples, int nSamples) {
////        const PoissonState<AGENT> &prior = windows.size()==0?startStatePrior:windows.back().analysis;
////        const ModelState<AGENT> &startState = windows.size()==0 ? startStatePrior.nextSample() : windows.back().realTrajectory.endState();
////        Trajectory<AGENT> realTrajectory(nTimesteps, startState);
//        windows.emplace_back(nTimesteps, std::move(analysisPMF), std::move(analysisSampler), pMakeObservation, pObserveIfPresent);
//        MCMCSolver<AGENT> solver(windows.back());
//        solver.solve(nBurnInSamples, nSamples);
//        analysis = std::move(solver.exactEndState);
//        analysisPMF = analysis.PMF();
//        analysisSampler = analysis.sampler();
//        return windows.back();
//    }
////
////
////    // Calculates the information gained by assimilation over nTimesteps
////    std::vector<double> calculateInformationGain() {
////        std::vector<double> informationGains(windows.size());
////        std::vector<PoissonState<AGENT>> referencePriors = priorWindows(10000);
////        for(int w=0; w<windows.size(); ++w) {
////            informationGains[w] = informationGain(
////                    windows[w].realTrajectory.endState(),
////                    referencePriors[w],
////                    windows[w].analysis);
////            Gnuplot gp;
////            Experiments::plotHeatMap(gp, referencePriors[w], windows[w].realTrajectory.endState());
////            Gnuplot gp2;
////            Experiments::plotHeatMap(gp2, windows[w].analysis, windows[w].realTrajectory.endState());
////        }
////        return informationGains;
////    }
////
////
//////    static PoissonState<AGENT> assimilateWindow(
//////            int nTimeseps,
//////            std::vector<Observation<AGENT>> observations,
//////            const PoissonState<AGENT> &startStatePrior,
////////            const Trajectory<AGENT> &initialTrajectory,
//////            int nSamples) {
//////
//////        PredPreyProblem<AGENT> abm(nTimeseps, observations, [&](const Trajectory<AGENT> &exactEndState) {
//////            return startStatePrior.extendedLogProb(exactEndState(0));
//////        });
//////        SimplexMCMC mcmc(abm, abm.logProbFunc());
//////        Trajectory<AGENT> initialTrajectory = generateInitialState(mcmc, startStatePrior);
//////        mcmc.setLPState(initialTrajectory);
//////        assert(mcmc.abmSanityChecks());
//////
//////        // test
//////        std::vector<int> occupationHistogram(10,0.0);
//////
//////        PoissonState<AGENT> finalState;
//////        for(int n=0; n<nSamples; ++n) {
//////            mcmc.nextSample();
//////            if(nSamples<1000 || n%1000 == 1) {
////////                std::cout << "Sample " << n << std::endl;
//////                assert(abm.isValidSolution(mcmc.X()));
////////            std::cout << "Sample " << n << " : " << glp::SparseVec(mcmc.X()) << std::endl;
//////            }
//////            const Trajectory<AGENT> &exactEndState = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
//////            ModelState<AGENT> endStateSample = exactEndState(abm.nTimesteps);
//////            occupationHistogram[endStateSample[AGENT(152)]] += 1;
//////            finalState += endStateSample;
//////        }
//////        debug(
//////                std::cout << "Occupation histogram " << occupationHistogram << std::endl;
//////                std::cout << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl;
//////                std::cout << "Feasible nextSample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl;
//////                std::cout << "Infeasible nextSample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
//////        );
//////        return finalState;
//////    }
////
////
////    // Calculates the posiion states at the end of each window, given no observations but
////    // including the trajectoryPrior over the start state. This is needed in order to calculate the
////    // information gain due to assimilation.
////    std::vector<PoissonState<AGENT>> priorWindows(int nSamples) {
////        int sumOfTimesteps = 0;
////        for(const AssimilationWindow<AGENT> &window: windows) {
////            sumOfTimesteps += window.realTrajectory.nTimesteps();
////        }
////        std::vector<PoissonState<AGENT>> endOfWindowStates(windows.size());
////        for(int s=0; s<nSamples; ++s) {
////            Trajectory<AGENT> trajectory(sumOfTimesteps, startStatePrior.nextSample());
////            int t = 0;
////            for(int w=0; w<windows.size(); ++w) {
////                t += windows[w].realTrajectory.nTimesteps();
////                endOfWindowStates[w] += trajectory(t);
////            }
////        }
////        return endOfWindowStates;
////    }
////
////
//////    static std::vector<PoissonState<AGENT>> posterior(
//////            int nWindows,
//////            int nTimestepsPerWindow,
//////            int nSamplesPerWindow,
//////            const PoissonState<AGENT> &startState,
//////            const std::vector<std::vector<Observation<AGENT>>> &observations) {
//////        ModelState<AGENT> realState = startState.nextSample();
//////        std::vector<PoissonState<AGENT>> endOfWindowStates(nWindows);
//////
//////        const PoissonState<AGENT> *windowStart = &startState;
//////        for (int window = 0; window < nWindows; ++window) {
//////            endOfWindowStates[window] = DataAssimilation<AGENT>::assimilateWindow(
//////                    nTimestepsPerWindow,
//////                    observations,
//////                    *windowStart,
//////                    nSamplesPerWindow
//////            );
//////            windowStart = &endOfWindowStates.back();
//////        }
//////        return endOfWindowStates;
//////    }
////
////
//////    static std::vector<std::vector<Observation<AGENT>>> observe(
//////            const Trajectory<AGENT> &realTrajectory,
//////            int nWindows,
//////            double pMakeObservation,
//////            double pObserveIfPresent) {
//////        std::vector<std::vector<Observation<AGENT>>> observations(nWindows);
//////
//////        const ModelState<AGENT> *windowStartState = &startState;
//////        for (int window = 0; window < nWindows; ++window) {
//////            auto obs = Observation<AGENT>::generateObservations(*windowStartState, nTimestepsPerWindow,
//////                                                                                 pMakeObservation, pObserveIfPresent);
//////            observations.push_back(obs);
//////            endOfWindowStates.push_back(exactEndState(nTimestepsPerWindow));
//////            windowStartState = &endOfWindowStates.back();
//////        }
//////        return std::pair(observations, endOfWindowStates);
//////    }
////
////
////    static double  informationGain(const ModelState<AGENT> &realState, const PoissonState<AGENT> &prior, const PoissonState<AGENT> &posterior) {
////        return (posterior.extendedLogProb(realState) - prior.extendedLogProb(realState))/log(2);
////    }
////
////
//};

#endif //GLPKTEST_DATAASSIMILATION_H
