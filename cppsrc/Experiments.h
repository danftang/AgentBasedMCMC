//
// Created by daniel on 06/05/2021.
//

#ifndef GLPKTEST_EXPERIMENTS_H
#define GLPKTEST_EXPERIMENTS_H


#include "gnuplot-iostream/gnuplot-iostream.h"
#include "StateTrajectory.h"
#include "agents/PredPreyAgent.h"
#include "PoissonState.h"
#include "BinomialDistribution.h"
#include "IntSampleStatistics.h"

class Experiments {
public:

//    static void PredPreyExpt();
    static void CatMouseExpt();
    static void RandomWalk();
//    static void GnuplotTest();

    static double nullPMF(const std::vector<double> &X) { return 0.0; }

//    static Gnuplot &plotHeatMap(Gnuplot &gp, const PoissonState<PredPreyAgent> &aggregateState, const ModelState<PredPreyAgent> &realState);
    static Gnuplot &plotAgents(Gnuplot &gp, const ModelState<PredPreyAgent> &state);

    static void PredPreyAssimilation();

    static std::vector<double> informationIncrease(int argc, char **argv);

    static std::vector<double>
    informationIncrease(int gridsize, int windowSize, int nWindows, double pPredator, double pPrey,
                        double pMakeObservation, double pObserveIfPreset, int nSamplesPerWindow, int nBurnInSamples);

//    Gnuplot &
//    plotHeatMap(Gnuplot &gp, const BinomialDistribution &aggregateState, const ModelState<PredPreyAgent> &realState);
    static void CatMouseAssimilation();

    static void BinomialAgentAssimilation();

    template<typename DOMAIN>
    static double informationGain(
            const DOMAIN &realState,
            const IntSampleStatistics<DOMAIN> &prior,
            const IntSampleStatistics<DOMAIN> &analysis) {
        for(int i=0; i < realState.size(); ++i) {
            std::cout << "Real occupancy = " << realState[i]
                      << " prior p = " << prior.P(i, realState[i])
                      << " posterior p = " << analysis.P(i, realState[i])
                      << " post / prior p = " << analysis.P(i, realState[i]) / prior.P(i, realState[i])
                      << std::endl;
        }
        return (analysis.logP(realState) - prior.logP(realState)) / log(2.0);

    }

    static void PredPreyExpt();
};


#endif //GLPKTEST_EXPERIMENTS_H
