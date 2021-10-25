#include <iostream>
#include <thread>
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

using glp::X;

//
//class MyClass {
//public:
//    int i;
//
//    MyClass(): i(0) { std::cout << "Default constructing" << std::endl;}
//    MyClass(int j): i(j) { std::cout << "Initialising" << std::endl; }
//    MyClass(const MyClass &other): i(other.i) { std::cout << "Copying" << std::endl;}
//    MyClass(MyClass &&other): i(other.i) { other.i = -1; std::cout << "Moving" << std::endl;}
//
//    int operator()() const { return i; }
//
//    operator std::function<int()>() const & { return [*this]() { return i; }; }
//    operator std::function<int()>() && { return [c = MyClass(std::move(*this))]() { return c.i; }; }
//
//    void myFunc(const MyClass &other) { std::cout << "myFunc const lValue ref" << std::endl;}
//    void myFunc(MyClass &&other) { std::cout << "myFunc rValue ref" << std::endl;}
//
//    template<typename T>
//    void myTFunc(T &other) { std::cout << "myTFunc lValue ref" << std::endl;}
//    template<typename T>
//    void myTFunc(T &&other) { std::cout << "myTFunc rValue ref" << std::endl;}
//};

//template<typename T>
//class MyClass {
//public:
//    T i;// = 1234;
//
//    MyClass(T p) { i = p; }
//
//    T x() { return i; }
//
//    int myFunc();
//
////    double logP(const std::vector<double> &X) { return i; }
////
////    int operator()() { return i; }
////    double operator()(double x) { return i+x; }
//
//
//};

//template<class... Types>
//void f(Types... inits) {
////    std::tuple<Types...> t(inits...);
////    std::index_sequence_for<Types...> indices;
//
//    ((std::cout << inits << std::endl),...);
//
//}
//
//template<class... Types, size_t... Indices>
//void expandTuple(std::tuple<Types...> &tuple, std::index_sequence<Indices...> indx) {
//    ((std::cout << std::get<Indices>(tuple) << std::endl),...);
//}
//
//template<class... Types>
//void expandTuple(std::tuple<Types...> &tuple) {
//    expandTuple(tuple,std::index_sequence_for<Types...>());
//}


int main(int argc, char *argv[]) {
//    Experiments::BinomialAgentAssimilation();

    using namespace dataflow;

    dataflow::Producer<int> p = []() { return 5; };
    dataflow::Consumer<int> c = [](int x) { std::cout << x << std::endl; return true; };

//    p >>= Take(5) >>= c;
//    p >>= Take(5) >>= Map<int(int)>([](int x) { return x+1;}) >>= c;
    p >>= Take{5} >>= Split {
        Take{2} >>= Map<int(int)>{[](int x) { return x+1;}} >>= c,
        switchAfter(
                3,
                Map<int(int)>([](int x) { return x+2;}) >>= c,
                Map<int(int)>([](int x) { return x+200;}) >>= c)
    };

//    p >>= Split {
//        Take(5) >>= Map([](int x) { return x+1; }) >>= c,
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

