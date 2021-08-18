//
// Created by daniel on 23/07/2021.
//

#ifndef GLPKTEST_UNITTESTS_H
#define GLPKTEST_UNITTESTS_H


#include <vector>
#include <boost/math/distributions/binomial.hpp>
//#include "ActFermionicDistribution.h"
//#include "StlStream.h"
//#include "ConvexPMF.h"
//#include "agents/BinomialAgent.h"
//#include "TrajectoryPriorPMF.h"
//#include "TrajectoryLikelihoodPMF.h"
//#include "TrajectoryPriorSampler.h"
//#include "DeltaPMF.h"
//#include "PoissonPMF.h"
//#include "BinomialDistribution.h"
//#include "AssimilationProblem.h"
//#include "RejectionSampler.h"
//#include "BinarySolutionSet.h"
#include "AssimilationWindow.h"
#include "agents/BinomialAgent.h"
#include "BinomialDistribution.h"
#include "BinarySolutionSet.h"
#include "RejectionSampler.h"

class UnitTests: public AssimilationWindow<BinomialAgent> {
public:
//    int nTimesteps;
//    BinomialDistribution startDist;
//    AgentStateObservation<BinomialAgent>    observation;
//    TrajectoryPriorPMF<BinomialAgent>       priorPMF;
//    TrajectoryPriorSampler<BinomialAgent>   priorSampler;
//    TrajectoryLikelihoodPMF<BinomialAgent>  likelihoodPMF;
//    ConvexPMF                               posterior;

//    TrajectoryPriorDistribution<BinomialAgent> window;
//    AssimilationProblem problem;


    UnitTests()
    :
            AssimilationWindow<BinomialAgent>(
                    2,
                    BinomialDistribution({
                                                 boost::math::binomial_distribution(1.0, 0.9),
                                                 boost::math::binomial_distribution(1.0, 0.1),
                                                 boost::math::binomial_distribution(1.0, 0.0)
                                         }),
                    AgentStateObservation<BinomialAgent>(
                            State<BinomialAgent>(1, 0),
                            1,
                            0.9
                    )
            )
    {
    }

//    static AssimilationProblem
//    initProblem(const TrajectoryPriorDistribution<BinomialAgent> &window, const BinomialDistribution &startDist, const AgentStateObservation<BinomialAgent> &observation) {
//        BinomialAgent::GRIDSIZE = window.nTimesteps + 1;
//        AssimilationProblem problem(
//                window.trajectoryPrior(startDist.PMF()),
//                window.priorSampler(startDist.sampler())
//        );
//        problem.addObservation(window.likelihood(observation));
//        return problem;
//    }


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


    void testValidTrajectorySet() {
        ConvexPMF pmf = priorPMF;
        ModelState<BinomialAgent> marginalisedFinalState;
        double marginalP = 0.0;
        for(const std::vector<double> &traj: BinarySolutionSet(pmf.convexSupport, pmf.nDimensions)) {
            double jointP = exp(pmf.logP(traj));
            marginalP += jointP;
//            std::cout << traj << " " << jointP << std::endl;
            ModelState<BinomialAgent> endState = Trajectory<BinomialAgent>(traj).endState();
            endState *= jointP;
            marginalisedFinalState += endState;
        }
        std::cout << "Marginal P = " << marginalP << std::endl;
        marginalisedFinalState *= 1.0/marginalP;
        std::cout << marginalisedFinalState << std::endl;
    }


    void testABMPrior() {
        ModelState<BinomialAgent> finalState;
        auto sampler = priorSampler;
        for(int s=0; s<100000; ++s) {
            Trajectory<BinomialAgent> sample = sampler();
//            std::cout << "Got sample: " << s << " " << sample << std::endl;
            finalState += sample.endState();
        }
        std::cout << "Prior samples: " <<  finalState << std::endl;
    }


    void testSimplexSampler() {
        ConvexPMF pmf = posterior; // window.trajectoryPrior(startDist.PMF());
//        std::cout << "Start state dist is:\n" << startDist.weights << std::endl;
//        std::cout << "Start state convexSupport is:\n" << startDist.PMF().convexSupport << std::endl;
        std::cout << "Simplex Sampler convexSupport:\n" << pmf.convexSupport << std::endl;

        SimplexMCMC sampler = SimplexMCMC(pmf, priorSampler());
//        std::cout << "kProbTokSim = " << sampler.kProbTokSim << std::endl;
//        std::cout << "kSimTokProb = " << sampler.kSimTokProb << std::endl;
        for(int burnIn=0; burnIn<1000; ++burnIn) {
            sampler.nextSample();
        }
        std::cout << "Sampler:\n" << sampler << std::endl;
        std::cout << "Initial solution: " << sampler.X() << std::endl;
        ModelState<BinomialAgent> mcmcFinalState;
        const int NSAMPLES = 1000000;
        for(int s=0; s<NSAMPLES; ++s) {
            Trajectory<BinomialAgent> sample(sampler.nextSample());
//            std::cout << "Sampler:\n" << sampler << std::endl;
//            std::cout << "Sample: " << sample << std::endl;
            assert(pmf.convexSupport.isValidSolution(sample));
            mcmcFinalState += sample.endState();
        }
        std::cout << "Feasible stats:\n" << sampler.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats:\n" << sampler.infeasibleStatistics << std::endl;
        std::cout << "Infeasible proportion = " << sampler.infeasibleStatistics.nSamples*1.0/(sampler.feasibleStatistics.nSamples + sampler.infeasibleStatistics.nSamples) << std::endl;
        mcmcFinalState *= 1.0/NSAMPLES;
        std::cout << mcmcFinalState << std::endl;
    }


    void testRejectionSampler() {
        std::cout << "Prior Support:\n" << priorPMF.convexSupport << std::endl;
        std::cout << "Likelihood Support:\n" << likelihoodPMF.convexSupport << std::endl;
        std::cout << "Convex support: \n" << posterior.convexSupport << std::endl;

        RejectionSampler rejectionSampler(priorSampler, likelihoodPMF);
        ModelState<BinomialAgent> finalState;
        const int NSAMPLES = 100000;
        for (int s = 0; s < NSAMPLES; ++s) {
            Trajectory<BinomialAgent> sample = rejectionSampler();
            finalState += sample.endState();
        }
        finalState *= 1.0/NSAMPLES;
        std::cout << "Mean state: " << finalState << std::endl;
    }


};


#endif //GLPKTEST_UNITTESTS_H
