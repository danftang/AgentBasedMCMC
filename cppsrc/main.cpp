#include <iostream>
#include <thread>
#include <future>
#include "glpkppinclude/glpkpp.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"
#include "UnitTests.h"
#include "TrajectoryPriorDistribution.h"
#include "ConvexPolyhedron.h"
#include "ABMConstraints.h"
#include "AssimilationProblem.h"
#include "diagnostics/Dataflow.h"
#include "diagnostics/MeanAndVariance.h"
#include "Plotter.h"

using glp::X;

int main(int argc, char *argv[]) {

//    using namespace dataflow;
//    auto m = [](int x) { return x*2; };
//    auto p = [n=0]() mutable { return ++n; };
//    auto c = [](const auto &x) { std::cout << x << std::endl; return true; };
//
//    p >>= Drop(10) >>= Take(1000) >>= Split {
//        Map(m) >>= c,
//        Map([](int x) { return x*x; }) >>= c
//    };

//    Experiments::DataflowDemo();
//    Experiments::PredPreyConvergence();

//    Experiments::BinomialAgentAssimilation();
//  Experiments::CatMouseSingleObservation();
//    Experiments::CatMouseAssimilation();
//    Experiments::CatMouseMultiObservation();
//    Experiments::PredPreySingleObservation();
//    Experiments::PredPreyAssimilation();

//    Experiments::FermionicIntegrality();

//    PoissonState<PredPreyAgent> startState;
//    TrajectoryPriorDistribution<PredPreyAgent> myPrior(startState, 8);

//    UnitTests tests;
//    tests.testPriorSampler();
//    tests.testRejectionSampler();
//    tests.testABMPrior();
//    tests.testSimplexSampler();
//    tests.testExactSolver();

//        std::cout << Experiments::informationIncrease(
//            8,
//            2,
//            1,
//            0.16,
//            0.32,
//            0.05,
//            0.9,
//            100000,
//            1000) << std::endl;


//    Experiments::RandomWalk();
//    Experiments::GnuplotTest();

    return 0;
}

