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
#include "ConstrainedFactorisedSampler.h"

template<int GRIDSIZE, int TIMESTEPS>
class FiguresForPaper {
public:
    typedef PredPreyTrajectory<GRIDSIZE,TIMESTEPS> trajectory_type;
    typedef typename trajectory_type::value_type value_type;

    static const std::string filenamePrefix;
    static const std::string problemFilename;
    static const std::string statFilename;


    static void generateStandardProblemFile(double kappa)  {
        constexpr double pPredator = PredPreyAgent<GRIDSIZE>::pPred;
        constexpr double pPrey = PredPreyAgent<GRIDSIZE>::pPrey;
        constexpr double pMakeObservation = 0.001;//0.05;
        constexpr double pObserveIfPresent = 1.0;
        std::ofstream probFile(problemFilename);
        boost::archive::binary_oarchive probArchive(probFile);

        Random::gen.seed(1234);
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem(pPredator, pPrey, pMakeObservation, pObserveIfPresent, kappa);
        std::cout << problem << std::endl;
        probArchive << problem;
    }


    static void generateStats(int nSamples) {
        constexpr int nThreads = 4;

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
        probArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

//        std::cout << "Basis vectors:" << std::endl;
//        for(const auto &basisVec : problem.basisObj.basisVectors) std::cout << basisVec << std::endl;

        auto posterior = problem.posterior();

        auto startTime = std::chrono::steady_clock::now();

        std::future<std::vector<ChainStats>> futureResults[nThreads];
        for(int thread = 0; thread < nThreads; ++thread) {
            futureResults[thread] = std::async(&startStatsThread,
                                               posterior,
                                               problem.basisVectors(),
                                               problem.randomInitialSolution(),
                                               problem.kappa,
                                               nSamples);
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
    }


    static std::vector<ChainStats> startStatsThread(
            const ConstrainedFactorisedDistribution<trajectory_type> &targetDistribution,
            const std::vector<SparseVec<value_type>> &basisVectors,
            trajectory_type initalSample,
            double kappa,
            int nSamples) {
        using namespace dataflow;

        const double maxLagProportion = 0.5;
        const int nLags = 200;
        const int nBurnIn = nSamples*0.2;
        ABM::kappa = kappa;

        ConstrainedFactorisedSampler sampler(targetDistribution, basisVectors, initalSample);

        std::valarray<std::valarray<double>> firstSynopsisSamples(nSamples/2);
        std::valarray<std::valarray<double>> lastSynopsisSamples(nSamples/2);
        ModelState<PredPreyAgent<GRIDSIZE>> firstMeanEndState;
        ModelState<PredPreyAgent<GRIDSIZE>> lastMeanEndState;

        sampler //>>= [](const trajectory_type &item) -> const trajectory_type & { std::cout << item << std::endl; return item; }
                >>= Drop(nBurnIn)
                >>= &PredPreyTrajectory<GRIDSIZE,TIMESTEPS>::endState
                >>= SwitchOnClose {
                        Split {
                                synopsis >>= save(firstSynopsisSamples),
                                Sum(nSamples/2,firstMeanEndState)
                        },
                        Split{
                                synopsis >>= save(lastSynopsisSamples),
                                Sum(nSamples/2,firstMeanEndState)
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


    static void singleThreadStats(int nSamples) {
        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
        probArchive >> problem;
        ABM::kappa = problem.kappa;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

        auto posterior = problem.posterior();

        auto startTime = std::chrono::steady_clock::now();

        std::vector<ChainStats> stats = startStatsThread(posterior,
                                                         problem.basisVectors(),
                                                         problem.randomInitialSolution(),
                                                         problem.kappa,
                                                         nSamples);

        auto endTime = std::chrono::steady_clock::now();
//        std::cout << stats << std::endl;
        std::cout << "Executed in " << endTime - startTime << std::endl;
    }


    static void sampleTiming(int nSamples) {
        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
        probArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

        std::cout << "Creating modelStateSampler..." << std::endl;
        ConstrainedFactorisedSampler sampler(problem.tableau, problem.posterior.factors, problem.posterior.perturbableFunctionFactory(), problem.prior.nextSample(false), problem.kappa);

        std::cout << "Burning in..." << std::endl;
        for(int burnIn=0; burnIn<nSamples/4; ++burnIn) sampler();

        std::cout << "Taking samples..." << std::endl;
        auto startTime = std::chrono::steady_clock::now();
        for(int s=0; s<nSamples; ++s) sampler();
        auto endTime = std::chrono::steady_clock::now();
        auto execTimeMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

        const trajectory_type &sample = sampler();
        std::cout << "Final state is " << std::vector(sample.begin(), sample.begin() + 20) << "..." << std::endl;
        std::cout << "Executed in " << execTimeMilliSeconds << "ms" << std::endl;
        std::cout << "Time per feasible sample " << execTimeMilliSeconds*1.0/nSamples << "ms" << std::endl;

    }


    static void plotProblemEndState() {
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem(problemFilename);
        Plotter plotter;
        plotter.plot(ModelState(problem.realTrajectory, problem.realTrajectory.nTimesteps()));
    }



    static void plotStats(bool waitForKeypressToExit = false) {
        std::ifstream statFile(statFilename);
        if(!statFile.good()) throw("Can't open stats probFile. Maybe you haven't run the analysis for this geometry yet.");
        boost::archive::binary_iarchive statArchive(statFile);
        MultiChainStats stats(0.0,"");
        statArchive >> stats;

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Odd because the stats probFile exists!");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
        probArchive >> problem;

        std::cout << problem << std::endl;

        // plot autocorrelations
        Plotter acPlotter;
        std::valarray<std::valarray<double>> autocorrelation = stats.autocorrelation();
        int nDimensions = autocorrelation[0].size();
        double xStride = stats.front().varioStride;
        acPlotter << "plot ";
        for(int d=1; d<=nDimensions; ++d) acPlotter << "'-' using (" << xStride <<  "*$0):" << d << " with lines title 'Statistic " << d << "', ";
        acPlotter << "0 with lines notitle\n";
        for(int d=1; d<=nDimensions; ++d) acPlotter.send1d(autocorrelation);

        // Print MCMC stats
        std::cout << stats;

        // Print scale reduction and effective samples
        std::valarray<double> neff = stats.effectiveSamples();
        std::valarray<double> ineff = (stats.nSamplesPerChain() * 1.0) / neff;
        double execTimePerSample = stats.execTimeMilliSeconds * 2.0 / (stats.front().stats.totalProposals() * stats.size());
        double execTimePerFeasibleSample = stats.execTimeMilliSeconds * 1.0 / (stats.nSamplesPerChain() * stats.nChains());
        std::cout << "Summary statistics for " << GRIDSIZE << " x " << TIMESTEPS << std::endl;
        std::cout << "Total exec time: " << stats.execTimeMilliSeconds << "ms" << std::endl;
        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
        std::cout << "Actual number of samples per chain: " << stats.nSamplesPerChain() << std::endl;
        std::cout << "Number of chains: " << stats.size() << std::endl;
        std::cout << "Effective number of samples (per chain): " << neff << std::endl;
        std::cout << "Sample inefficiency factor: " << ineff << std::endl << std::endl;
        std::cout << "Execution time per sample (all threads): " << execTimePerSample << "ms" << std::endl;
        std::cout << "Execution time per feasible sample (all threads): " << execTimePerFeasibleSample << "ms" << std::endl;
        std::cout << "Execution time per (worst case) effective sample: " << stats.execTimeMilliSeconds/(neff.min()*stats.nChains()) << "ms" << std::endl;
        // plot end state

        Plotter endStatePlotter;
        endStatePlotter.plot(problem.realTrajectory.endState(), stats.meanEndState(),"");

        if(waitForKeypressToExit) {
            std::cout << "Press Enter to exit" << std::endl;
            std::cin.get();
        }

    }


    static std::valarray<double> synopsis(const ModelState<PredPreyAgent<GRIDSIZE>> &endState) {
        std::valarray<double> synopsis(floor(log2(GRIDSIZE)) -1);
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
            origin += partitionSize;
        }
        return synopsis;
    }

};

template<int GRIDSIZE, int TIMESTEPS>
const std::string FiguresForPaper<GRIDSIZE,TIMESTEPS>::filenamePrefix = "../data/PredPrey" + std::to_string(GRIDSIZE) + "-" + std::to_string(TIMESTEPS);

template<int GRIDSIZE, int TIMESTEPS>
const std::string FiguresForPaper<GRIDSIZE,TIMESTEPS>::problemFilename = FiguresForPaper<GRIDSIZE,TIMESTEPS>::filenamePrefix + ".prob";

template<int GRIDSIZE, int TIMESTEPS>
const std::string FiguresForPaper<GRIDSIZE,TIMESTEPS>::statFilename = FiguresForPaper<GRIDSIZE,TIMESTEPS>::filenamePrefix + ".stat";


#endif //GLPKTEST_FIGURESFORPAPER_H
