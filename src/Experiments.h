//
// Created by daniel on 06/05/2021.
//

#ifndef GLPKTEST_EXPERIMENTS_H
#define GLPKTEST_EXPERIMENTS_H

#include <functional>
#include "Trajectory.h"
#include "RejectionSampler.h"
#include "Forecast.h"
#include "Likelihood.h"
#include "diagnostics/AgentDataflow.h"

using namespace dataflow;

class Experiments {
public:
    static void BinomialAgentSingleObservation();
    static void CatMouseSingleObservation();
    static void PredPreySingleObservation();

    static void CatMouseAssimilation();

public:

//    template<class STARTSTATESAMPLER, class AGENT>
//    static void doRejectionSample(
//            STARTSTATESAMPLER &startstateSampler,
//            TrajectoryForecastDistribution<AGENT> &forecast,
//            ObservationLikelihood<AGENT> &likelihood,
//            int nSamples);

    template<class AGENT>
    static void doSingleObservationExperiment(
            int nTimesteps,
            int nBurnin,
            int nSamples,
            int nRejectionSamples,
            double kappa,
            double alpha,
            const StartStateDistribution<AGENT> &startState,
            const AgentStateObservation<AGENT> &observation);


////    static void PredPreySingleObservation();
//    static void CatMouseMultiObservation();
//    static void RandomWalk();
////    static void GnuplotTest();
//
//    static double nullPMF(const std::vector<double> &X) { return 0.0; }
//
////    static Gnuplot &plotHeatMap(Gnuplot &gp, const PoissonState<PredPreyAgent> &aggregateState, const ModelState<PredPreyAgent> &realState);
////    static Gnuplot &plotAgents(Gnuplot &gp, const ModelState<PredPreyAgent> &state);
//
//    static void PredPreyAssimilation();
//
////    static std::vector<double> informationIncrease(int argc, char **argv);
//
//    template<int GRIDSIZE>
//    static std::vector<double>
//    informationIncrease(int windowSize, int nWindows, double pPredator, double pPrey,
//                        double pMakeObservation, double pObserveIfPreset, int nSamplesPerWindow, int nBurnInSamples);
//
////    Gnuplot &
////    plotHeatMap(Gnuplot &gp, const BinomialDistribution &aggregateState, const ModelState<PredPreyAgent> &realState);
//
//
//    template<typename DOMAIN>
//    static double informationGain(
//            const DOMAIN &realState,
//            const IntSampleStatistics<DOMAIN> &prior,
//            const IntSampleStatistics<DOMAIN> &analysis) {
////        for(int i=0; i < realState.size(); ++i) {
////            std::cout << "Real occupancy = " << realState[i]
////                      << " prior p = " << prior.P(i, realState[i])
////                      << " posterior p = " << analysis.P(i, realState[i])
////                      << " post / prior p = " << analysis.P(i, realState[i]) / prior.P(i, realState[i])
////                      << std::endl;
////        }
//        return (analysis.logP(realState) - prior.logP(realState)) / log(2.0);
//
//    }
//
//
//    static void CatMousePrior();
//
//    static void FermionicIntegrality();
//
//    template<int GRIDSIZE>
//    static std::valarray<double> Synopsis(const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory);
//
//    template<int GRIDSIZE>
//    static MultiChainStats PredPreyConvergenceThread(const ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> &posterior,
//                              Trajectory<PredPreyAgent<GRIDSIZE>> startState);
//
//    static void PredPreyConvergence();
//
//    static void DataflowDemo();
//
//    static void animatedPredPreyDemo();
//
//    static void minimalBasis();
};


#endif //GLPKTEST_EXPERIMENTS_H
