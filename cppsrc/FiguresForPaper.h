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
    static void generateProblemFile();

    template<int GRIDSIZE>
    static void generateStats(const char *problemFile, const char *outputFile) {
        constexpr int nThreads = 4;

        std::ifstream file(problemFile);
        boost::archive::binary_iarchive boostArchive(file);
        PredPreyProblem<GRIDSIZE> problem;
        boostArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;


        std::future<MultiChainStats> futureResults[nThreads];
        for(int thread = 0; thread < nThreads; ++thread) {
            futureResults[thread] = std::async(&startStatsThread<GRIDSIZE>, problem.posterior(), problem.priorSampler()());
        }

        MultiChainStats multiChainStats;
        multiChainStats.reserve(2*nThreads);
        for(int thread=0; thread<nThreads; ++thread) {
            futureResults[thread].wait();
            multiChainStats += futureResults[thread].get();
        }

        std::ofstream outFile(outputFile);
        boost::archive::binary_oarchive boostout(outFile);

        boostout << multiChainStats;
        std::cout << multiChainStats;
    }

    template<int GRIDSIZE>
    static MultiChainStats startStatsThread(ConvexPMF<Trajectory<PredPreyAgent<GRIDSIZE>>> posterior, Trajectory<PredPreyAgent<GRIDSIZE>> startState) {
        using namespace dataflow;
        constexpr int nSamples = 50000; // must be an even number
        assert((nSamples&1) == 0);
        const double maxLagProportion = 0.4;
        const int nLags = 80;
        constexpr int nBurnIn = nSamples*0.5;
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
        std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*100.0/sampler.simplex.feasibleStatistics.nSamples << "%" << std::endl;

        firstMeanEndState /= nSamples/2;
        lastMeanEndState /= nSamples/2;

        MultiChainStats stats;
        stats.reserve(2);
        stats.emplace_back(std::move(firstSynopsisSamples), nLags, maxLagProportion, std::move(firstMeanEndState), std::move(firstNextSample[0]));
        stats.emplace_back(std::move(lastSynopsisSamples), nLags, maxLagProportion, std::move(lastMeanEndState), std::move(lastNextSample[0]));
//    std::cout << stats << std::endl;
        return std::move(stats);
    }

    template<int GRIDSIZE>
    static void plotStats(const char *statsFile, const char *problemFile) {
        MultiChainStats stats;
        std::ifstream inFile(statsFile);
        boost::archive::binary_iarchive boostin(inFile);
        boostin >> stats;

        std::ifstream file(problemFile);
        boost::archive::binary_iarchive boostArchive(file);
        PredPreyProblem<GRIDSIZE> problem;
        boostArchive >> problem;

        std::cout << stats;

        Plotter gp;
        std::valarray<std::valarray<double>> autocorrelation = stats.autocorrelation();
        gp.heatmap(autocorrelation, 0.5, autocorrelation[0].size()-0.5, 0.5*stats.front().varioStride, (autocorrelation.size()+0.5)*stats.front().varioStride);

        Plotter gp2;
        gp2 << "plot '-' using 0:1 with lines, '-' using 0:2 with lines, 0 with lines\n";
        gp2.send1d(autocorrelation);
        gp2.send1d(autocorrelation);

        std::valarray<double> neff = stats.effectiveSamples();
        std::valarray<double> ineff = (stats.nSamples()*1.0)/neff;

        std::cout << stats << std::endl;
        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
        std::cout << "Actual number of samples per chain: " << stats.nSamples() << std::endl;
        std::cout << "Effective number of samples: " << neff << std::endl;
        std::cout << "Sample inefficiency factor: " << ineff << std::endl;

        Plotter gp3;
        gp3.plot(problem.realTrajectory.endState(), stats.meanEndState());


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

};


#endif //GLPKTEST_FIGURESFORPAPER_H
