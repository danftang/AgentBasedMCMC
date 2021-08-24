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

    static double informationGain(const std::vector<double> &realEndState, const BinomialDistribution &prior,
                           const BinomialDistribution &analysis);
};


#endif //GLPKTEST_EXPERIMENTS_H
