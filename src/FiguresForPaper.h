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
#include "FactorisedDistributionSampler.h"

template<int GRIDSIZE, int TIMESTEPS>
class FiguresForPaper {
public:
    typedef PredPreyAgent<GRIDSIZE> agent_type;
    typedef PredPreyTrajectory<GRIDSIZE,TIMESTEPS> trajectory_type;
    typedef BernoulliStartState<agent_type> start_state_type;
    typedef typename trajectory_type::value_type value_type;

    static const std::string filenamePrefix;
    static const std::string problemFilename;
    static const std::string statFilename;


    static void generateStandardPredPreyPosteriorFile(double kappa)  {
        constexpr double pPredator = PredPreyAgent<GRIDSIZE>::pPred;
        constexpr double pPrey = PredPreyAgent<GRIDSIZE>::pPrey;
        constexpr double pMakeObservation = 0.001;//0.05;
        constexpr double pObserveIfPresent = 1.0;
        std::ofstream probFile(problemFilename);
        boost::archive::binary_oarchive probArchive(probFile);

        Random::gen.seed(123456);
        start_state_type startState([](agent_type agent) {
            return agent.type()==PredPreyAgent<GRIDSIZE>::PREDATOR?pPredator:pPrey;
        }, kappa);
//        ABMPrior prior = makeABMPrior<trajectory_type>(startState);
//        ABMLikelihood likelihood(prior.nextSample(), pMakeObservation, pObserveIfPresent, kappa);
        ABMPosterior posterior = makeABMPosterior<trajectory_type>(startState, pMakeObservation, pObserveIfPresent, kappa);

        std::cout << posterior << std::endl;
        probArchive << posterior;
    }


    static void generateStats(int nSamples, int nThreads = 4) {

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        ABMPosterior<trajectory_type, start_state_type> posterior;
        probArchive >> posterior;

        std::cout << "Loaded ABM posterior:" << std::endl;
        std::cout << posterior;
//        std::cout << "Basis vectors:" << std::endl;
//        for(const auto &basisVec : problem.basisObj.basisVectors) std::cout << basisVec << std::endl;

        FactorisedDistributionSampler sampler(posterior);
        //
//        std::future<std::pair<ChainStats,ChainStats>> futureResults[nThreads];
//        for(int thread = 0; thread < nThreads; ++thread) {
//            futureResults[thread] = std::async(&startStatsThread,
//                                               posterior,
//                                               problem.randomInitialSolution(),
//                                               nSamples);
//        }
//
//        MultiChainStats multiChainStats;
//        multiChainStats.reserve(2*nThreads);
//        for(int thread=0; thread<nThreads; ++thread) {
//            futureResults[thread].wait();
//            multiChainStats.add(futureResults[thread].get());
//        }

        MultiChainStats multiChainStats = MultiChainStats::analyseConvergence(sampler, nSamples, nThreads);

        std::ofstream statFile(statFilename);
        if(!statFile.good()) throw("Can't open stats probFile to save results.");
        boost::archive::binary_oarchive statArchive(statFile);
        statArchive << multiChainStats;
    }


//    static std::pair<ChainStats,ChainStats> startStatsThread(
//            const ConstrainedFactorisedDistribution<trajectory_type> &targetDistribution,
//            trajectory_type initalSample,
//            int nSamples) {
//        using namespace dataflow;
//
//        const double maxLagProportion = 0.5;
//        const int nLags = 200;
//        const int nBurnIn = nSamples*0.2;
//
//        FactorisedDistributionSampler sampler(targetDistribution, initalSample);
//
//        std::valarray<std::valarray<double>> firstSynopsisSamples(nSamples/2);
//        std::valarray<std::valarray<double>> lastSynopsisSamples(nSamples/2);
//        ModelState<PredPreyAgent<GRIDSIZE>> meanEndState;
//
//        for(int s = 0; s<nBurnIn; ++s) sampler();
//        for (int s = 0; s < nSamples; ++s) {
//            auto endState = sampler().endState();
//            meanEndState += endState;
//            if (s < nSamples / 2)
//                firstSynopsisSamples[s - nBurnIn] = synopsis(endState);
//            else
//                lastSynopsisSamples[s - nBurnIn] = synopsis(endState);
//        }
//
//        std::cout << "Stats =\n" << sampler.stats << std::endl;
//
//        std::valarray<double> meanEndState1 = meanEndState / nSamples;
//
//        std::pair<ChainStats, ChainStats> stats(
//                ChainStats(
//                        firstSynopsisSamples,
//                        nLags,
//                        maxLagProportion,
//                        sampler.stats),
//                ChainStats(
//                        lastSynopsisSamples,
//                        nLags,
//                        maxLagProportion,
//                        sampler.stats)
//        );
//        return stats;
//    }


//    static void singleThreadStats(int nSamples) {
//        std::ifstream probFile(problemFilename);
//        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
//        boost::archive::binary_iarchive probArchive(probFile);
//        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
//        probArchive >> problem;
//        ABM::kappa = problem.kappa;
//
//        std::cout << "Loaded problem" << std::endl;
//        std::cout << problem;
//
//        auto posterior = problem.posterior();
//
//        auto startTime = std::chrono::steady_clock::now();
//
//        std::vector<ChainStats> stats = startStatsThread(posterior,
//                                                         problem.basisVectors(),
//                                                         problem.randomInitialSolution(),
//                                                         problem.kappa,
//                                                         nSamples);
//
//        auto endTime = std::chrono::steady_clock::now();
////        std::cout << stats << std::endl;
//        std::cout << "Executed in " << endTime - startTime << std::endl;
//    }


//    static void sampleTiming(int nSamples) {
//        std::ifstream probFile(problemFilename);
//        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
//        boost::archive::binary_iarchive probArchive(probFile);
//        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
//        probArchive >> problem;
//
//        std::cout << "Loaded problem" << std::endl;
//        std::cout << problem;
//
//        std::cout << "Creating modelStateSampler..." << std::endl;
//        FactorisedDistributionSampler sampler(problem.tableau, problem.posterior.factors, problem.posterior.perturbableFunctionFactory(), problem.prior.nextSample(false), problem.kappa);
//
//        std::cout << "Burning in..." << std::endl;
//        for(int burnIn=0; burnIn<nSamples/4; ++burnIn) sampler();
//
//        std::cout << "Taking samples..." << std::endl;
//        auto startTime = std::chrono::steady_clock::now();
//        for(int s=0; s<nSamples; ++s) sampler();
//        auto endTime = std::chrono::steady_clock::now();
//        auto execTimeMilliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
//
//        const trajectory_type &sample = sampler();
//        std::cout << "Final state is " << std::vector(sample.begin(), sample.begin() + 20) << "..." << std::endl;
//        std::cout << "Executed in " << execTimeMilliSeconds << "ms" << std::endl;
//        std::cout << "Time per feasible sample " << execTimeMilliSeconds*1.0/nSamples << "ms" << std::endl;
//
//    }


//    static void plotProblemEndState() {
//        ABMPosterior<trajectory_type,start_state_type> posterior(problemFilename);
//        Plotter plotter;
//        plotter.plot(ModelState(posterior.likelihood.realTrajectory.value());
//    }



    static void plotStats(bool waitForKeypressToExit = false) {
        std::ifstream statFile(statFilename);
        if(!statFile.good()) throw("Can't open stats probFile. Maybe you haven't run the analysis for this geometry yet.");
        boost::archive::binary_iarchive statArchive(statFile);
        MultiChainStats stats;
        statArchive >> stats;

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Odd because the stats probFile exists!");
        boost::archive::binary_iarchive probArchive(probFile);
        ABMPosterior<trajectory_type,start_state_type> posterior;
        probArchive >> posterior;

        std::cout << posterior << std::endl;

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
        auto execTimePerSample = stats.execTime * 2.0 / (stats.front().samplerStats.totalProposals() * stats.size());
        auto execTimePerFeasibleSample = stats.execTime * 1.0 / (stats.nSamplesPerChain() * stats.nChains());
        std::cout << "Summary statistics for " << GRIDSIZE << " x " << TIMESTEPS << std::endl;
        std::cout << "Total exec time: " << stats.execTime << std::endl;
        std::cout << "Potential scale reduction: " << stats.potentialScaleReduction() << std::endl;
        std::cout << "Actual number of samples per chain: " << stats.nSamplesPerChain() << std::endl;
        std::cout << "Number of chains: " << stats.size() << std::endl;
        std::cout << "Effective number of samples (per chain): " << neff << std::endl;
        std::cout << "Sample inefficiency factor: " << ineff << std::endl << std::endl;
        std::cout << "Execution time per sample (all threads): " << execTimePerSample << std::endl;
        std::cout << "Execution time per feasible sample (all threads): " << execTimePerFeasibleSample << std::endl;
        std::cout << "Execution time per (worst case) effective sample: " << stats.execTime/(neff.min()*stats.nChains()) << std::endl;
        // plot end state

        Plotter endStatePlotter;
        endStatePlotter.plot(posterior.likelihood.realTrajectory.value().endState(), stats.meanEndState(),"");

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
