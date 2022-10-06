//
// Created by daniel on 03/10/22.
//

#ifndef ABMCMC_CONVERGENCEANALYSIS_H
#define ABMCMC_CONVERGENCEANALYSIS_H

#include <vector>
#include <iostream>
#include <valarray>
#include <future>
#include "ChainStats.h"
#include "../FactorisedDistributionSampler.h"
#include "../ModelState.h"
#include "MultiChainStats.h"

template<class DOMAIN>
class ConvergenceAnalysis {
public:
    typedef DOMAIN trajectory_type;
    typedef typename DOMAIN::agent_type agent_type;

    static void generateStats(int nSamples) {
        constexpr int nThreads = 4;

        std::ifstream probFile(problemFilename);
        if(!probFile.good()) throw("Can't open problem probFile for this geometry. Run generateStandardProblemFile first.");
        boost::archive::binary_iarchive probArchive(probFile);
        PredPreyProblem<GRIDSIZE,TIMESTEPS> problem;
        probArchive >> problem;

        std::cout << "Loaded problem" << std::endl;
        std::cout << problem;

        auto posterior = problem.posterior();

        auto startTime = std::chrono::steady_clock::now();

        std::future<std::tuple<ChainStats,ChainStats,ModelState<PredPreyAgent<GRIDSIZE>>>> futureResults[nThreads];
        for(int thread = 0; thread < nThreads; ++thread) {
            futureResults[thread] = std::async(&startStatsThread,
                                               posterior,
                                               posterior.basis.randomInitialSolution(),
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
            const FactorisedDistributionOnBasis<trajectory_type> &targetDistribution,
            trajectory_type initalSample,
            int nSamples) {

        const double maxLagProportion = 0.5;
        const int nLags = 200;
        const int nBurnIn = nSamples*0.2;

        FactorisedDistributionSampler sampler(targetDistribution, initalSample);

        std::valarray<std::valarray<double>> firstSynopsisSamples(nSamples/2);
        std::valarray<std::valarray<double>> lastSynopsisSamples(nSamples/2);
        ModelState<agent_type> firstMeanEndState;
        ModelState<agent_type> lastMeanEndState;

        for(int s = 0; s<nBurnIn; ++s) sampler();
        for (int s = 0; s < nSamples; ++s) {
            auto endState = sampler().endState();
            if (s < nSamples / 2) {
                firstSynopsisSamples[s - nBurnIn] = synopsis(endState);
                firstMeanEndState += endState;
            } else {
                lastSynopsisSamples[s - nBurnIn] = synopsis(endState);
                lastMeanEndState += endState;
            }
        }

        std::cout << "Stats =\n" << sampler.stats << std::endl;

        std::valarray<double> meanEndState1 = firstMeanEndState / nSamples;

        std::vector<ChainStats> stats;
        stats.reserve(2);
        stats.emplace_back(
                firstSynopsisSamples,
                nLags,
                maxLagProportion,
                sampler.stats);
        stats.emplace_back(
                lastSynopsisSamples,
                nLags,
                maxLagProportion,
                sampler.stats);
        //    std::cout << stats << std::endl;
        return stats;
    }

};


#endif //ABMCMC_CONVERGENCEANALYSIS_H
