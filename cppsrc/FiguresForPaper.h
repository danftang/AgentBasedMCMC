//
// Created by daniel on 09/11/2021.
//

#ifndef GLPKTEST_FIGURESFORPAPER_H
#define GLPKTEST_FIGURESFORPAPER_H

#include <future>
#include <thread>
#include "PredPreyProblem.h"
#include "diagnostics/MultiChainStats.h"
#include "diagnostics/Dataflow.h"
#include "Plotter.h"
#include "MCMCSampler.h"

class FiguresForPaper {
public:

    template<int GRIDSIZE>
    static void generateStatsAndPlot(int nTimesteps) {
        generateStats<GRIDSIZE>(nTimesteps);
        plotStats<GRIDSIZE>(nTimesteps);
    }

    static void generateAllProblemFiles() {
        for (int nTimesteps = 2; nTimesteps <= 16; nTimesteps *= 2) {
            generateStandardProblemFile<8>(nTimesteps);
            generateStandardProblemFile<16>(nTimesteps);
            generateStandardProblemFile<32>(nTimesteps);
        }
    }

    template<int GRIDSIZE>
    static void generateStandardProblemFile(int nTimesteps)  {
        constexpr double pPredator = 0.05; // (0.08 + 10.24/(GRIDSIZE*GRIDSIZE))/3.0;
        constexpr double pPrey = 0.05; // 2.0*pPredator;
        constexpr double pMakeObservation = 0.02;
        constexpr double pObserveIfPresent = 0.9;
        std::string probFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
        std::ofstream probFile(probFilename);
        boost::archive::binary_oarchive probArchive(probFile);

        Random::gen.seed(1234);
        PredPreyProblem<GRIDSIZE> problem(nTimesteps, pPredator, pPrey, pMakeObservation, pObserveIfPresent);
        std::cout << problem << std::endl;
        probArchive << problem;
    }


    template<int GRIDSIZE>
    static void generateStats(int nTimesteps) {
        std::string problemFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
        std::string statFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".stat";
        constexpr int nThreads = 4;

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE> problem;
        probArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

        auto startTime = std::chrono::steady_clock::now();

        std::future<MultiChainStats> futureResults[nThreads];
        for(int thread = 0; thread < nThreads; ++thread) {
            futureResults[thread] = std::async(&startStatsThread<GRIDSIZE>, problem.posterior(), problem.priorSampler()());
        }


        MultiChainStats multiChainStats(problemFilename);
        multiChainStats.reserve(2*nThreads);
        for(int thread=0; thread<nThreads; ++thread) {
            futureResults[thread].wait();
            multiChainStats += futureResults[thread].get();
        }

        auto endTime = std::chrono::steady_clock::now();

        multiChainStats.execTimeMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        std::ofstream statFile(statFilename);
        if(!statFile.good()) throw("Can't open stats probFile to save results.");
        boost::archive::binary_oarchive statArchive(statFile);
        statArchive << multiChainStats;
//        std::cout << std::endl << "Saved stats" << std::endl;
//        std::cout << multiChainStats;
    }

    template<int GRIDSIZE>
    static MultiChainStats startStatsThread(ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> posterior, Trajectory<PredPreyAgent<GRIDSIZE>> startState) {
        using namespace dataflow;
        constexpr int nSamples = 4000000; // must be an even number
        assert((nSamples&1) == 0);
        const double maxLagProportion = 0.5;
        const int nLags = 100;
        constexpr int nBurnIn = 200000;//nSamples*0.25;
        auto trajectoryToEnergy = [](const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory) { return -trajectory.logProb(); };
        auto trajectoryToEndState = [](const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory) { return trajectory.endState(); };
        MCMCSampler sampler(posterior, startState);

        std::valarray<std::valarray<double>> firstSynopsisSamples(nSamples/2);
        std::valarray<std::valarray<double>> lastSynopsisSamples(nSamples/2);
        ModelState<PredPreyAgent<GRIDSIZE>> firstMeanEndState;
        ModelState<PredPreyAgent<GRIDSIZE>> lastMeanEndState;
        std::valarray<Trajectory<PredPreyAgent<GRIDSIZE>>> firstNextSample(Trajectory<PredPreyAgent<GRIDSIZE>>(0),1);
        std::valarray<Trajectory<PredPreyAgent<GRIDSIZE>>> lastNextSample(Trajectory<PredPreyAgent<GRIDSIZE>>(0),1);

        sampler >>= Drop(nBurnIn)
                >>= SwitchOnClose {
                        Split {
                                Map{ synopsis<GRIDSIZE> } >>= save(firstSynopsisSamples),
                                Take(nSamples/2) >>= Map{ trajectoryToEndState } >>= Sum(firstMeanEndState)
                        },
                        save(firstNextSample),
                        Split{
                                Map{ synopsis<GRIDSIZE> } >>= save(lastSynopsisSamples),
                                Take(nSamples/2) >>= Map{ trajectoryToEndState } >>= Sum(firstMeanEndState)
                        },
                        save(lastNextSample)
                };

        std::cout << "Feasible stats =\n" << sampler.simplex.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats =\n" << sampler.simplex.infeasibleStatistics << std::endl;
        std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/(sampler.simplex.feasibleStatistics.nSamples + sampler.simplex.infeasibleStatistics.nSamples) << "%" << std::endl << std::endl;

        firstMeanEndState /= nSamples/2;
        lastMeanEndState /= nSamples/2;

        MultiChainStats stats;
        stats.reserve(2);
        stats.emplace_back(
                std::move(firstSynopsisSamples),
                nLags,
                maxLagProportion,
                std::move(firstMeanEndState),
                std::move(firstNextSample[0]),
                sampler.simplex.feasibleStatistics,
                sampler.simplex.infeasibleStatistics);
        stats.emplace_back(
                std::move(lastSynopsisSamples),
                nLags,
                maxLagProportion,
                std::move(lastMeanEndState),
                std::move(lastNextSample[0]),
                sampler.simplex.feasibleStatistics,
                sampler.simplex.infeasibleStatistics);
        //    std::cout << stats << std::endl;
        return std::move(stats);
    }

    template<int GRIDSIZE>
    static void plotStats(int nTimesteps) {
        std::string statFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".stat";
        std::ifstream statFile(statFilename);
        if(!statFile.good()) throw("Can't open stats probFile. Maybe you haven't run the analysis for this geometry yet.");
        boost::archive::binary_iarchive statArchive(statFile);
        MultiChainStats stats;
        statArchive >> stats;

        std::string probFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
        std::ifstream probFile(probFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Odd because the stats probFile exists!");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE> problem;
        probArchive >> problem;

        // plot autocorrelations
        Plotter acPlotter;
        std::valarray<std::valarray<double>> autocorrelation = stats.autocorrelation();
        int nDimensions = autocorrelation[0].size();
        double xStride = stats.front().varioStride;
        acPlotter << "set title 'Autocorrelations " << GRIDSIZE << " x " << nTimesteps << "'\n";
        acPlotter << "plot ";
        for(int d=1; d<=nDimensions; ++d) acPlotter << "'-' using (" << xStride <<  "*$0):" << d << " with lines title 'synopsis " << d << "', ";
        acPlotter << "0 with lines notitle\n";
        for(int d=1; d<=nDimensions; ++d) acPlotter.send1d(autocorrelation);

        // Print MCMC stats
        std::cout << stats;
//        for(const ChainStats &chain: stats) {
//            std::cout << "Feasible MCMC stats:" << std::endl;
//            std::cout << chain.feasibleStats << std::endl;
//            std::cout << "Infeasible MCMC stats:" << std::endl;
//            std::cout << chain.infeasibleStats << std::endl;
//            std::cout << "Infeasible proportion = "
//            << chain.infeasibleStats.nSamples*100.0/(chain.feasibleStats.nSamples + chain.infeasibleStats.nSamples)
//            << "%" << std::endl << std::endl;
//        }

        // Print scale reduction and effective samples
        std::valarray<double> neff = stats.effectiveSamples();
        std::valarray<double> ineff = (stats.nSamples()*1.0)/neff;
        double execTimePerSample = stats.execTimeMilliSeconds * 2.0 / (stats.front().feasibleStats.nSamples*stats.size());
        std::cout << "Summary statistics for " << GRIDSIZE << " x " << nTimesteps << std::endl;
        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
        std::cout << "Actual number of samples per chain: " << stats.nSamples() << std::endl;
        std::cout << "Effective number of samples: " << neff << std::endl;
        std::cout << "Sample inefficiency factor: " << ineff << std::endl << std::endl;
        std::cout << "Execution time per sample: " << execTimePerSample << "ms" << std::endl;
        std::cout << "Mean execution time per effective sample: " << ineff.sum()*execTimePerSample/ineff.size() << "ms" << std::endl;
        // plot end state

        Plotter().plot(problem.realTrajectory.endState(), stats.meanEndState(),"End state " + std::to_string(GRIDSIZE) + " x " + std::to_string(nTimesteps));
    }

    template<int GRIDSIZE>
    static std::valarray<double> synopsis(const Trajectory<PredPreyAgent<GRIDSIZE>> &trajectory) {
        std::valarray<double> synopsis(floor(log2(GRIDSIZE)) -1);
        ModelState<PredPreyAgent<GRIDSIZE>> endState = trajectory.endState();
        int origin = 0;
        int varid = 0;
        for(int partitionSize = GRIDSIZE / 2; partitionSize > 1; partitionSize /=2) {
            double predOccupation = 0.0;
            double preyOccupation = 0.0;
            for(int x=0; x < partitionSize; ++x) {
                for(int y=0; y < partitionSize; ++y) {
                    predOccupation += endState[PredPreyAgent<GRIDSIZE>(origin + x,origin + y,PredPreyAgent<GRIDSIZE>::PREDATOR)];
                    preyOccupation += endState[PredPreyAgent<GRIDSIZE>(origin + x,origin + y,PredPreyAgent<GRIDSIZE>::PREY)];
                }
            }
            assert(varid < synopsis.size());
            synopsis[varid++] = predOccupation + preyOccupation;
//        synopsis[varid++] = preyOccupation;
            origin += partitionSize;
        }
        return synopsis;
    }

    static std::string filenamePrefix(int gridsize, int timesteps) {
        return "PredPrey" + std::to_string(gridsize) + "-" + std::to_string(timesteps);
    }

};


#endif //GLPKTEST_FIGURESFORPAPER_H
