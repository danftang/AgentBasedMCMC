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
//#include "TrajectoryLikelihoodPMF.h"
//#include "TrajectorySampler.h"
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
#include "ExactSolver.h"
#include "BernoulliModelState.h"
#include "MCMCSampler.h"

class UnitTests: public AssimilationWindow<BinomialAgent> {
public:
//    int nTimesteps;
//    BinomialDistribution startDist;
//    AgentStateObservation<BinomialAgent>    observation;
//    TrajectoryPMF<BinomialAgent>       priorPMF;
//    TrajectorySampler<BinomialAgent>   priorSampler;
//    TrajectoryLikelihoodPMF<BinomialAgent>  likelihoodPMF;
//    ConvexPMF                               posterior;

//    TrajectoryPriorDistribution<BinomialAgent> window;
//    AssimilationProblem problem;


    UnitTests()
            :
            AssimilationWindow<BinomialAgent>(
                    initProblem(
                            2,
                            BernoulliModelState<BinomialAgent>({0.9, 0.1, 0.0}),
                            AgentStateObservation<BinomialAgent>(
                                    State<BinomialAgent>(1, 0),
                                    1,
                                    0.9
                            )
                    )
            ) {
    }

    static AssimilationWindow<BinomialAgent>
    initProblem(int nTimesteps, const BernoulliModelState<BinomialAgent> &startDist, const AgentStateObservation<BinomialAgent> &observation) {
        BinomialAgent::GRIDSIZE = nTimesteps + 1;
        return AssimilationWindow<BinomialAgent>(
                2,
                startDist,
                observation);
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
        ConvexPMF<std::vector<double>> myPMF([](const std::vector<double> &X) {
            return log(X[1] + X[2] + X[3]);
        },4);

        myPMF.convexSupport.push_back(0 <= 1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 2);
        myPMF.convexSupport.push_back(0 <= 1.0*X(1) <= 1);
        myPMF.convexSupport.push_back(0 <= 1.0*X(2) <= 1);
        myPMF.convexSupport.push_back(0 <= 1.0*X(3) <= 1);
        std::cout << myPMF.convexSupport << std::endl;

        std::vector<int> sampleCounts(8,0);
        MCMCSampler sampler(myPMF);
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


    void testExactSolver() {
        ExactSolver<BinomialAgent> solver(posterior);
        std::cout << "Exact exactEndState = " << solver.exactEndState << std::endl;
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

        MCMCSampler sampler(pmf, priorSampler());
//        std::cout << "kProbTokSim = " << sampler.kProbTokSim << std::endl;
//        std::cout << "kSimTokProb = " << sampler.kSimTokProb << std::endl;
        for(int burnIn=0; burnIn<1000; ++burnIn) {
            sampler.nextSample();
        }
        std::cout << "Sampler:\n" << sampler.simplex << std::endl;
        std::cout << "Initial exactEndState: " << sampler.simplex.X() << std::endl;
        ModelState<BinomialAgent> mcmcFinalState;
        const int NSAMPLES = 1000000;
        for(int s=0; s<NSAMPLES; ++s) {
            Trajectory<BinomialAgent> sample(sampler.nextSample());
//            std::cout << "Sampler:\n" << sampler << std::endl;
//            std::cout << "Sample: " << sample << std::endl;
            assert(pmf.convexSupport.isValidSolution(sample));
            mcmcFinalState += sample.endState();
        }
        std::cout << "Feasible stats:\n" << sampler.simplex.feasibleStatistics << std::endl;
        std::cout << "Infeasible stats:\n" << sampler.simplex.infeasibleStatistics << std::endl;
        std::cout << "Infeasible proportion = " << sampler.simplex.infeasibleStatistics.nSamples*1.0/(sampler.simplex.feasibleStatistics.nSamples + sampler.simplex.infeasibleStatistics.nSamples) << std::endl;
        mcmcFinalState *= 1.0/NSAMPLES;
        std::cout << mcmcFinalState << std::endl;
    }


    void testRejectionSampler() {
        std::cout << "Prior Support:\n" << priorPMF.convexSupport << std::endl;
        std::cout << "Likelihood Support:\n" << likelihoodPMF.convexSupport << std::endl;
        std::cout << "Convex support: \n" << posterior.convexSupport << std::endl;

        RejectionSampler rejectionSampler(priorSampler, likelihoodPMF);
        ModelState<BinomialAgent> finalState;
        const int NSAMPLES = 200000;
        for (int s = 0; s < NSAMPLES; ++s) {
            Trajectory<BinomialAgent> sample = rejectionSampler();
            finalState += sample.endState();
        }
        finalState *= 1.0/NSAMPLES;
        std::cout << "Mean state: " << finalState << std::endl;

        ExactSolver<BinomialAgent> exactSolver(posterior);
        std::cout << "Exact solution: " << exactSolver.exactEndState << std::endl;

    }

    void testPriorSampler() {
        ModelState<BinomialAgent> finalState;
        const int NSAMPLES = 100000;
        for (int s = 0; s < NSAMPLES; ++s) {
            Trajectory<BinomialAgent> sample = priorSampler();
//            std::cout << "Sample = " << sample << std::endl;
//            std::cout << "Sample start state = " << sample(0) << std::endl;
//            std::cout << "Sample end state = " << sample.endState() << std::endl;
            finalState += sample.endState();
        }
        finalState *= 1.0/NSAMPLES;
        std::cout << "Sampler mean: " << finalState << std::endl;

        ExactSolver<BinomialAgent> exactSolver(priorPMF);
        std::cout << "Exact solution: " << exactSolver.exactEndState << std::endl;

    }

    static void testConstraintImplication() {
        CatMouseAgent leftMouse(CatMouseAgent::MOUSE, CatMouseAgent::LEFT);
        Event<CatMouseAgent> leftMouseMove(0, leftMouse, CatMouseAgent::MOVE);
        Event<CatMouseAgent> leftMouseStay(0, leftMouse, CatMouseAgent::STAYPUT);

        std::vector<glp::Constraint> constraints;
        std::vector<glp::Constraint> y = leftMouse.constraints(0, CatMouseAgent::MOVE);
        ABMConstraints<CatMouseAgent>::push_xImpliesY(constraints, leftMouseMove, y[0]);

        std::cout << leftMouseMove << " implies " << y[0] << std::endl;
        std::cout << constraints << std::endl << std::endl;

        constraints.clear();
        y.clear();
        y = leftMouse.constraints(0, CatMouseAgent::STAYPUT);
        ABMConstraints<CatMouseAgent>::push_xImpliesY(constraints, leftMouseStay, y[0]);

        std::cout << leftMouseStay << " implies " << y[0] << std::endl;
        std::cout << constraints << std::endl;

    }


};


#endif //GLPKTEST_UNITTESTS_H
