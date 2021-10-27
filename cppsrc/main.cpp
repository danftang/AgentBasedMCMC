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

using glp::X;

template<typename C, typename R, typename A>
void unwrap(R(C::*f)(A) const) {
}

auto doThread(int nBurnIn, int nSamples) {
    using namespace dataflow;
    Producer<int> sampler = [n = 0]() mutable { return ++n; };
    std::vector<double> energy;
    Consumer<int> c = [](int x) { std::cout << x << std::endl; return true; };
    auto trajectoryToEnergy = [](int i) { return -i*0.1; };
    auto synopsis = [](int i) { return std::vector<double>{ 1.0*i }; };
    MeanAndVariance  meanVariances;

    sampler >>= Drop(nBurnIn) >>= Take(nSamples) >>= Split {
            Map { trajectoryToEnergy } >>= pushBack(energy),
            Map { synopsis } >>= meanVariances.consumer()
    };
    return std::pair(energy,meanVariances);
}

int main(int argc, char *argv[]) {
//    Experiments::BinomialAgentAssimilation();

//    using namespace dataflow;
//    int nBurnIn = 100;
//    int nSamples = 200;
//    Producer<int> sampler = [n = 0]() mutable { return ++n; };
//    std::vector<double> energy;
//    Consumer<int> c = [](int x) { std::cout << x << std::endl; return true; };
//    auto trajectoryToEnergy = [](int i) { return -i*0.1; };
//    auto synopsis = [](int i) { return std::vector<double>{ 1.0*i }; };
//    MeanAndVariance  meanVariances;
//
//    sampler >>= Drop(nBurnIn) >>= Take(nSamples) >>= Split {
//            Map { trajectoryToEnergy } >>= pushBack(energy),
//            Map { synopsis } >>= meanVariances.consumer()
//    };
//
//    std::cout << energy << std::endl;
//    std::cout << meanVariances.mean() << std::endl;
//    std::cout << meanVariances.sampleVariance() << std::endl;

    auto threadResult = std::async(&doThread,100,200);

    std::cout << threadResult.get();

//    p >>= Take{5} >>= Split {
//        Take{2} >>= add1 >>= c,
//        SwitchAfter {
//            3,
//            Drop(1) >>= add1 >>= add1 >>= c,
//            add1 >>= add1 >>= add1 >>= c
//        }
//    };



//    unwrap(&decltype(f)::operator());
//    myFunc(std::function(f));

//    std::cout << d << std::endl;
//    myFunc(1);
//    myFunc(MyClass());
//    myFunc(f);

//    decltype(myFunc(std::function(std::declval<MyClass>()))) d = 1.234;
//    std::cout << d << std::endl;

//    p >>= Split {
//        Take(5) >>= Map<int(int)>([](int x) { return x+1; }) >>= c,
//        Take(4) >>= c,
//        Take(3) >>= c
//    };

//    dataflow::map([](int x) -> int { return x+1; });

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

