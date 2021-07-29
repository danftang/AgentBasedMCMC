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
#include "ABMPrior.h"
#include "DeltaPMF.h"

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
        });

        myPMF.convexSupport.addConstraint(0 <= 1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 2);
        myPMF.convexSupport.addConstraint(0 <= 1.0*X(1) <= 1);
        myPMF.convexSupport.addConstraint(0 <= 1.0*X(2) <= 1);
        myPMF.convexSupport.addConstraint(0 <= 1.0*X(3) <= 1);
        std::cout << myPMF.convexSupport << std::endl;

        std::vector<int> sampleCounts(8,0);
        for(int s=0; s<100000; ++s) {
            std::vector<double> sample = myPMF.nextSample();
            int vertexId = sample[1] + sample[2]*2 + sample[3]*4;
            ++sampleCounts[vertexId];
//            std::cout << sample[3] << sample[2] << sample[1] << std::endl;
        }

        for(int id=0; id <8; ++id) {
            std::cout << id << " " << sampleCounts[id] << std::endl;
        }
    }

    static void testABMPrior() {
        const int nTimesteps = 2;
        BinomialAgent::GRIDSIZE = nTimesteps+1;
        ModelState<BinomialAgent> startState;
        startState[0] = 1;
        DeltaPMF startPMF(startState);

        ABMPrior<BinomialAgent,DeltaPMF> prior(startPMF, nTimesteps);

        std::cout << prior.convexSupport << std::endl;

        ModelState<BinomialAgent> finalState;
        for(int s=0; s<10000; ++s) {
            Trajectory<BinomialAgent> sample = prior.nextSample();
            finalState += sample.endState();
        }
        std::cout << finalState << std::endl;

        ModelState<BinomialAgent> mcmcFinalState;
        ConvexPMF &convexPmf = prior;
        for(int burnIn=0; burnIn<100; ++burnIn) {
            convexPmf.nextSample();
        }
        std::cout << "Sampler:\n" << *convexPmf.sampler << std::endl;
        for(int s=0; s<10000; ++s) {
            Trajectory<BinomialAgent> sample = convexPmf.nextSample();
//            std::cout << sample << std::endl;
            mcmcFinalState += sample.endState();
        }
        std::cout << "Feasible stats:\n" << convexPmf.sampler->feasibleStatistics << std::endl;
        std::cout << "Infeasible stats:\n" << convexPmf.sampler->infeasibleStatistics << std::endl;
        std::cout << mcmcFinalState << std::endl;

        boost::math::binomial binom = boost::math::binomial(nTimesteps, BinomialAgent::pMove);
        for(int i=0; i<=nTimesteps; ++i) {
            std::cout << i << " -> " << boost::math::pdf(binom, i) << std::endl;
        }
    }
};


#endif //GLPKTEST_UNITTESTS_H
