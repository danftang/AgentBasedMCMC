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

//template<typename C, typename R, typename A>
//void unwrap(R(C::*f)(A) const) {
//}

template<typename... T> decltype(auto) myFunc(T&&...t) {
    return std::forward_as_tuple<T&&...>(std::forward<T>(t)...);
}

template<typename... T>
class MyClass {
public:
        std::tuple<T...> tup;
    MyClass(std::tuple<T...> &&tuple): tup(std::move(tuple)) {}
    MyClass(const std::tuple<T...> &tuple): tup(tuple) {}
    MyClass(T...vals): tup(std::move(vals)...) {}

//    template<typename... R>
//    MyClass(R&&...r): MyClass(
//            std::forward_as_tuple<R&&...>(std::forward<R>(r)...)
//                    ) {
//    }
};

int main(int argc, char *argv[]) {


//    Experiments::DataflowDemo();
    Experiments::PredPreyConvergence();

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

