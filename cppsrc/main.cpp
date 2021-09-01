#include <iostream>
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


int main(int argc, char *argv[]) {


//   Experiments::BinomialAgentAssimilation();
//    Experiments::CatMouseExpt();
//        Experiments::CatMouseAssimilation();
//    Experiments::PredPreyExpt();
    Experiments::PredPreyAssimilation();


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

