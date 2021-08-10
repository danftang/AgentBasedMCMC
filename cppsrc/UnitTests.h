//
// Created by daniel on 23/07/2021.
//

#ifndef GLPKTEST_UNITTESTS_H
#define GLPKTEST_UNITTESTS_H


#include <vector>
#include <boost/math/distributions/binomial.hpp>
#include "ActFermionicDistribution.h"
#include "StlStream.h"
#include "ConvexPMF.h"
#include "agents/BinomialAgent.h"
#include "TrajectoryDistribution.h"
#include "DeltaPMF.h"
#include "PoissonPMF.h"
#include "BinomialPMF.h"
#include "AssimilationProblem.h"
#include "RejectionSampler.h"
#include "ValidTrajectorySet.h"

class UnitTests {
public:

    static void testAll() {
        testActFermionicDistribution();
    }


    static void testActFermionicDistribution() {
        std::vector<double> p = {0.8, 0.15, 0.05};
        ActFermionicDistribution myDist(p);

        std::vector<int> actCounts(3,0);
        for(int s = 0; s<1000000; ++s) {
            std::vector<bool> acts = myDist.sampleUnordered(2);
            if(acts[0]) ++actCounts[0];
            if(acts[1]) ++actCounts[1];
            if(acts[2]) ++actCounts[2];
        }
        std::cout << actCounts << std::endl;

        std::vector<int> actCounts2(3,0);
        std::discrete_distribution<int> dist = std::discrete_distribution<int>(p.begin(), p.end());
        for(int s = 0; s<1000000; ) {
            int act1 = dist(Random::gen);
            int act2 = dist(Random::gen);
            if(act1 != act2) {
                ++actCounts2[act1];
                ++actCounts2[act2];
                ++s;
            }
        }
        std::cout << actCounts2 << std::endl;
        assert(abs(actCounts[0]-actCounts2[0]) < 1000);
        assert(abs(actCounts[1]-actCounts2[1]) < 1000);
        assert(abs(actCounts[2]-actCounts2[2]) < 1000);
    }


    static void testConvexPMF() {
        using glp::X;
        ConvexPMF myPMF([](const std::vector<double> &X) {
            return log(X[1] + X[2] + X[3]);
        },4);

        myPMF.convexSupport.push_back(0 <= 1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 2);
        myPMF.convexSupport.push_back(0 <= 1.0*X(1) <= 1);
        myPMF.convexSupport.push_back(0 <= 1.0*X(2) <= 1);
        myPMF.convexSupport.push_back(0 <= 1.0*X(3) <= 1);
        std::cout << myPMF.convexSupport << std::endl;

        std::vector<int> sampleCounts(8,0);
        SimplexMCMC sampler(myPMF);
        for(int s=0; s<100000; ++s) {
            const std::vector<double> sample = sampler.nextSample();
            int vertexId = sample[1] + sample[2]*2 + sample[3]*4;
            ++sampleCounts[vertexId];
//            std::cout << sample[3] << sample[2] << sample[1] << std::endl;
        }

        for(int id=0; id <8; ++id) {
            std::cout << id << " " << sampleCounts[id] << std::endl;
        }
    }


    static void testValidTrajectorySet() {
        const int nTimesteps = 2;
        BinomialAgent::GRIDSIZE = nTimesteps+1;
        BinomialPMF startDist({0.75, 0.75, 0.75}, 1);
        TrajectoryDistribution<BinomialAgent> window(nTimesteps);
        AssimilationProblem problem(
                window.prior(startDist.PMF()),
                window.priorSampler(startDist.sampler()));
        for(const std::vector<double> &traj: ValidTrajectorySet(problem)) {
            std::cout << traj << std::endl;
        }
    }

    static void testABMPrior() {
        const int nTimesteps = 2;
        BinomialAgent::GRIDSIZE = nTimesteps+1;

        BinomialPMF startDist({0.75, 0.75, 0.75}, 1);

        std::cout << "Start state convexSupport is:\n" << startDist.PMF().convexSupport << std::endl;


        TrajectoryDistribution<BinomialAgent> window(nTimesteps);

        AssimilationProblem problem(
                window.prior(startDist.PMF()),
                window.priorSampler(startDist.sampler()));

        std::cout << "Prior convexSupport:\n" << problem.priorPMF.convexSupport << std::endl;

        ModelState<BinomialAgent> finalState;
        for(int s=0; s<100000; ++s) {
            Trajectory<BinomialAgent> sample = problem.priorSampler();
            finalState += sample.endState();
        }
        std::cout << "Prior samples: " <<  finalState << std::endl;

        SimplexMCMC sampler = problem.simplexSampler();
        for(int burnIn=0; burnIn<1000; ++burnIn) {
            sampler.nextSample();
        }
        std::cout << "Sampler:\n" << sampler << std::endl;
        ModelState<BinomialAgent> mcmcFinalState;
        for(int s=0; s<1000000; ++s) {
            Trajectory<BinomialAgent> sample(sampler.nextSample());
//            std::cout << sample << std::endl;
            mcmcFinalState += sample.endState();
        }
        std::cout << "Feasible stats:\n" << sampler.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats:\n" << sampler.infeasibleStatistics << std::endl;
        std::cout << mcmcFinalState << std::endl;

//        boost::math::binomial binom = boost::math::binomial(nTimesteps, BinomialAgent::pMove);
//        for(int i=0; i<=nTimesteps; ++i) {
//            std::cout << i << " -> " << boost::math::pdf(binom, i) << std::endl;
//        }
    }


    static void testRejectionSampler() {
        const int nTimesteps = 2;
        BinomialAgent::GRIDSIZE = nTimesteps + 1;

        BinomialPMF startDist({0.4, 0.01, 0.01}, 2);

        std::cout << "Start state convexSupport is:\n" << startDist.PMF().convexSupport << std::endl;

        TrajectoryDistribution<BinomialAgent> window(nTimesteps);

        AgentStateObservation<BinomialAgent> observation(State<BinomialAgent>(1,0), 1, 1.0);

        AssimilationProblem problem(
                window.prior(startDist.PMF()),
                window.priorSampler(startDist.sampler()));

//        problem.addObservation(window.likelihood(observation));

        std::cout << "Prior Support:\n" << problem.priorPMF.convexSupport << std::endl;
        std::cout << "Likelihood Support:\n" << problem.likelihoodPMF.convexSupport << std::endl;
        std::cout << "Convex support: \n" << problem.PMF().convexSupport << std::endl;


        RejectionSampler rejectionSampler(problem);
        ModelState<BinomialAgent> finalState;
        for (int s = 0; s < 100000; ++s) {
            Trajectory<BinomialAgent> sample = rejectionSampler();
            finalState += sample.endState();
        }
        std::cout << "Rejection samples: " << finalState << std::endl;

        ModelState<BinomialAgent> mcmcFinalState;
        SimplexMCMC mcmcSampler = problem.simplexSampler();
        for (int burnIn = 0; burnIn < 1000; ++burnIn) {
            mcmcSampler.nextSample();
        }
        std::cout << "Sampler:\n" << mcmcSampler << std::endl;
        for (int s = 0; s < 100000; ++s) {
            Trajectory<BinomialAgent> sample(mcmcSampler.nextSample());
            //            std::cout << sample << std::endl;
            mcmcFinalState += sample.endState();
        }
        std::cout << "Feasible stats:\n" << mcmcSampler.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats:\n" << mcmcSampler.infeasibleStatistics << std::endl;
        std::cout << mcmcFinalState << std::endl;

    }

    static void testABMPosterior() {
        const int nTimesteps = 2;
        BinomialAgent::GRIDSIZE = nTimesteps + 1;

        BinomialPMF startDist({0.4, 0.01, 0.01}, 2);

        std::cout << "Start state convexSupport is:\n" << startDist.PMF().convexSupport << std::endl;

        TrajectoryDistribution<BinomialAgent> window(nTimesteps);

        AgentStateObservation<BinomialAgent> observation(State<BinomialAgent>(1,0), 1, 1.0);

        AssimilationProblem problem(
                window.prior(startDist.PMF()),
                window.priorSampler(startDist.sampler()));

        problem.addObservation(window.likelihood(observation));

        std::cout << "Prior convexSupport:\n" << problem.priorPMF.convexSupport << std::endl;

        ModelState<BinomialAgent> finalState;
        for (int s = 0; s < 10000; ++s) {
            Trajectory<BinomialAgent> sample = problem.priorSampler();
            finalState += sample.endState();
        }
        std::cout << "Prior samples: " << finalState << std::endl;

        ModelState<BinomialAgent> mcmcFinalState;
        SimplexMCMC sampler = problem.simplexSampler();
        for (int burnIn = 0; burnIn < 100; ++burnIn) {
            sampler.nextSample();
        }
        std::cout << "Sampler:\n" << sampler << std::endl;
        for (int s = 0; s < 10000; ++s) {
            Trajectory<BinomialAgent> sample(sampler.nextSample());
            //            std::cout << sample << std::endl;
            mcmcFinalState += sample.endState();
        }
        std::cout << "Feasible stats:\n" << sampler.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats:\n" << sampler.infeasibleStatistics << std::endl;
        std::cout << mcmcFinalState << std::endl;

        //        boost::math::binomial binom = boost::math::binomial(nTimesteps, BinomialAgent::pMove);
        //        for(int i=0; i<=nTimesteps; ++i) {
        //            std::cout << i << " -> " << boost::math::pdf(binom, i) << std::endl;
        //        }
    }
};


#endif //GLPKTEST_UNITTESTS_H
