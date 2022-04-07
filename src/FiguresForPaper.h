//
// Created by daniel on 09/11/2021.
//

#ifndef GLPKTEST_FIGURESFORPAPER_H
#define GLPKTEST_FIGURESFORPAPER_H

#include <future>
#include <thread>
#include "PredPreyProblem.h"
#include "diagnostics/MultiChainStats.h"
#include "diagnostics/AgentDataflow.h"
#include "Plotter.h"
#include "SparseBasisSampler.h"

class FiguresForPaper {
public:

    template<int GRIDSIZE>
    static void generateStatsAndPlot(int nTimesteps, int nSamples) {
        generateStats<GRIDSIZE>(nTimesteps, nSamples);
        plotStats<GRIDSIZE>(nTimesteps);
    }

//    static void generateAllProblemFiles() {
//        for (int nTimesteps = 2; nTimesteps <= 16; nTimesteps *= 2) {
//            generateStandardProblemFile<8>(nTimesteps);
//            generateStandardProblemFile<16>(nTimesteps);
//            generateStandardProblemFile<32>(nTimesteps);
//        }
//    }

    template<int GRIDSIZE>
    static void generateStandardProblemFile(int nTimesteps, double kappa)  {
        constexpr double pPredator = 1.0-PredPreyAgent<GRIDSIZE>::pNoPred;
        constexpr double pPrey = 1.0-PredPreyAgent<GRIDSIZE>::pNoPrey;
        constexpr double pMakeObservation = 0.05;
        constexpr double pObserveIfPresent = 1.0;//0.9;
        std::string probFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
        std::ofstream probFile(probFilename);
        boost::archive::binary_oarchive probArchive(probFile);

        Random::gen.seed(1234);
        PredPreyProblem<GRIDSIZE> problem(nTimesteps, pPredator, pPrey, pMakeObservation, pObserveIfPresent, kappa);
        std::cout << problem << std::endl;
        probArchive << problem;
    }


    template<int GRIDSIZE>
    static void generateStats(int nTimesteps, int nSamples) {
        constexpr int nThreads = 4;

        std::string problemFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
        std::string statFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".stat";

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE> problem;
        probArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

        TableauNormMinimiser<ABM::occupation_type> basisTableau(problem.posterior.constraints);
//        std::cout << "Tableau = \n" << basisTableau << std::endl;

        auto startTime = std::chrono::steady_clock::now();

        std::future<std::vector<ChainStats>> futureResults[nThreads];
        for(int thread = 0; thread < nThreads; ++thread) {
            futureResults[thread] = std::async(&startStatsThread<GRIDSIZE>, problem, basisTableau, problem.prior.nextSample(), nSamples);
        }


        MultiChainStats multiChainStats(problem.kappa, problemFilename);
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
////        std::cout << std::endl << "Saved stats" << std::endl;
////        std::cout << multiChainStats;
    }

    template<int GRIDSIZE>
    static std::vector<ChainStats> startStatsThread(const PredPreyProblem<GRIDSIZE> &problem, const TableauNormMinimiser<ABM::occupation_type> &tableau, std::vector<ABM::occupation_type> initialState, int nSamples) {
        using namespace dataflow;

        const double maxLagProportion = 0.5;
        const int nLags = 100;
        const int nBurnIn = nSamples*0.1;
        int nTimesteps = problem.nTimesteps();

        SparseBasisSampler sampler(tableau, problem.posterior.factors, problem.posterior.perturbableFunctionFactory(), initialState, problem.kappa);

        std::valarray<std::valarray<double>> firstSynopsisSamples(nSamples/2);
        std::valarray<std::valarray<double>> lastSynopsisSamples(nSamples/2);
        ModelState<PredPreyAgent<GRIDSIZE>> firstMeanEndState;
        ModelState<PredPreyAgent<GRIDSIZE>> lastMeanEndState;

        sampler >>= Drop(nBurnIn)
                >>= TrajectoryToModelState<PredPreyAgent<GRIDSIZE>>(nTimesteps, nTimesteps)
                >>= SwitchOnClose {
                        Split {
                                Map{ synopsis<GRIDSIZE> } >>= save(firstSynopsisSamples),
                                Take(nSamples/2) >>= Sum<ModelState<PredPreyAgent<GRIDSIZE>>>(firstMeanEndState)
                        },
                        Split{
                                Map{ synopsis<GRIDSIZE> } >>= save(lastSynopsisSamples),
                                Take(nSamples/2) >>= Sum<ModelState<PredPreyAgent<GRIDSIZE>>>(firstMeanEndState)
                        },
                };

        std::cout << "Stats =\n" << sampler.stats << std::endl;

        std::valarray<double> meanEndState1 = firstMeanEndState / (nSamples/2);
        std::valarray<double> meanEndState2 = lastMeanEndState / (nSamples/2);

        std::vector<ChainStats> stats;
        stats.reserve(2);
        stats.emplace_back(
                std::move(firstSynopsisSamples),
                nLags,
                maxLagProportion,
                std::move(meanEndState1),
                sampler.stats);
        stats.emplace_back(
                std::move(lastSynopsisSamples),
                nLags,
                maxLagProportion,
                std::move(meanEndState2),
                sampler.stats);
        //    std::cout << stats << std::endl;
        return std::move(stats);
    }


    template<int GRIDSIZE>
    static void singleThreadStats(int nTimesteps, int nSamples) {

        std::string problemFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
        std::string statFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".stat";

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE> problem;
        probArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

        TableauNormMinimiser<ABM::occupation_type> basisTableau(problem.posterior.constraints);
//        std::cout << "Tableau = \n" << basisTableau << std::endl;

        auto startTime = std::chrono::steady_clock::now();

        std::vector<ChainStats> stats = startStatsThread<GRIDSIZE>(problem, basisTableau, problem.prior.nextSample(), nSamples);

        auto endTime = std::chrono::steady_clock::now();
        auto execTimeMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        std::cout << "Executed in " << execTimeMilliSeconds << "ms" << std::endl;
//        MultiChainStats multiChainStats(problem.kappa, problem.alpha, problemFilename);
//        multiChainStats.reserve(2);
//        multiChainStats += stats;
//
//        multiChainStats.execTimeMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
//
////        std::ofstream statFile(statFilename);
////        if(!statFile.good()) throw("Can't open stats probFile to save results.");
////        boost::archive::binary_oarchive statArchive(statFile);
////        statArchive << multiChainStats;
//
////        std::cout << std::endl << "Saved stats" << std::endl;
//        std::cout << multiChainStats;
    }


    template<int GRIDSIZE>
    static void plotStats(int nTimesteps) {
        std::string statFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".stat";
        std::ifstream statFile(statFilename);
        if(!statFile.good()) throw("Can't open stats probFile. Maybe you haven't run the analysis for this geometry yet.");
        boost::archive::binary_iarchive statArchive(statFile);
        MultiChainStats stats(0.0,"");
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
        double execTimePerSample = stats.execTimeMilliSeconds * 2.0 / (stats.front().stats.nSamples()*stats.size());
        std::cout << "Summary statistics for " << GRIDSIZE << " x " << nTimesteps << std::endl;
        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
        std::cout << "Actual number of samples per chain: " << stats.nSamples() << std::endl;
        std::cout << "Effective number of samples: " << neff << std::endl;
        std::cout << "Sample inefficiency factor: " << ineff << std::endl << std::endl;
        std::cout << "Execution time per sample: " << execTimePerSample << "ms" << std::endl;
        std::cout << "Mean execution time per effective sample: " << ineff.sum()*execTimePerSample/ineff.size() << "ms" << std::endl;
        // plot end state

        Plotter endStatePlotter;
        endStatePlotter.plot(problem.realTrajectory.endState(), stats.meanEndState(),"End state " + std::to_string(GRIDSIZE) + " x " + std::to_string(nTimesteps));

//        std::cout << "Press Enter to exit" << std::endl;
//        std::cin.get();

        //        while(acPlotter.is_open() || endStatePlotter.is_open());
    }

//    template<int GRIDSIZE>
//    static void plotStatsToFile(int nTimesteps) {
//        std::string statFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".stat";
//        std::ifstream statFile(statFilename);
//        if(!statFile.good()) throw("Can't open stats probFile. Maybe you haven't run the analysis for this geometry yet.");
//        boost::archive::binary_iarchive statArchive(statFile);
//        MultiChainStats stats;
//        statArchive >> stats;
//
//        std::string probFilename = filenamePrefix(GRIDSIZE, nTimesteps) + ".prob";
//        std::ifstream probFile(probFilename);
//        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Odd because the stats probFile exists!");
//        boost::archive::binary_iarchive probArchive(probFile);
//        PredPreyProblem<GRIDSIZE> problem;
//        probArchive >> problem;
//
//        // plot autocorrelations
//        Plotter acPlotter;
//        std::valarray<std::valarray<double>> autocorrelation = stats.autocorrelation();
//        int nDimensions = autocorrelation[0].size();
//        double xStride = stats.front().varioStride;
//        acPlotter << "set title 'Autocorrelations " << GRIDSIZE << " x " << nTimesteps << "'\n";
//        acPlotter << "plot ";
//        for(int d=1; d<=nDimensions; ++d) acPlotter << "'-' using (" << xStride <<  "*$0):" << d << " with lines title 'synopsis " << d << "', ";
//        acPlotter << "0 with lines notitle\n";
//        for(int d=1; d<=nDimensions; ++d) acPlotter.send1d(autocorrelation);
//
//        // Print MCMC stats
//        std::cout << stats;
//
//        // Print scale reduction and effective samples
//        std::valarray<double> neff = stats.effectiveSamples();
//        std::valarray<double> ineff = (stats.nSamples()*1.0)/neff;
//        double execTimePerSample = stats.execTimeMilliSeconds * 2.0 / (stats.front().feasibleStats.nSamples*stats.size());
//        std::cout << "Summary statistics for " << GRIDSIZE << " x " << nTimesteps << std::endl;
//        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
//        std::cout << "Actual number of samples per chain: " << stats.nSamples() << std::endl;
//        std::cout << "Effective number of samples: " << neff << std::endl;
//        std::cout << "Sample inefficiency factor: " << ineff << std::endl << std::endl;
//        std::cout << "Execution time per sample: " << execTimePerSample << "ms" << std::endl;
//        std::cout << "Mean execution time per effective sample: " << ineff.sum()*execTimePerSample/ineff.size() << "ms" << std::endl;
//        // plot end state
//
//        Plotter().plot(problem.realTrajectory.endState(), stats.meanEndState(),"End state " + std::to_string(GRIDSIZE) + " x " + std::to_string(nTimesteps));
//    }



    template<int GRIDSIZE>
    static std::valarray<double> synopsis(const ModelState<PredPreyAgent<GRIDSIZE>> &endState) {
        std::valarray<double> synopsis(floor(log2(GRIDSIZE)) -1);
//        ModelState<PredPreyAgent<GRIDSIZE>> endState(trajectory, trajectory.nTimesteps(), trajectory.nTimesteps()); //.endState();
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
        return "../data/PredPrey" + std::to_string(gridsize) + "-" + std::to_string(timesteps);
    }

};


#endif //GLPKTEST_FIGURESFORPAPER_H
