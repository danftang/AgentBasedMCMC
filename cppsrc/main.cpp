#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"
#include "UnitTests.h"
#include "TrajectoryDistribution.h"
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

//class MyClass {
//public:
//    int i = 1234;
//
//    double logP(const std::vector<double> &X) { return i; }
//
//    int operator()() { return i; }
//    double operator()(double x) { return i+x; }
//
//
//};

int main(int argc, char *argv[]) {

//    PoissonState<PredPreyAgent> startState;
//    TrajectoryDistribution<PredPreyAgent> myPrior(startState, 8);

    UnitTests tests;
//    tests.testRejectionSampler();
//    tests.testABMPrior();
    tests.testSimplexSampler();
    tests.testValidTrajectorySet();
    // TODO: Work out why adding BinomialPMF constraints breaks simplex sampling

    //    std::cout << Experiments::informationIncrease(
//            8,
//            2,
//            1,
//            0.16,
//            0.32,
//            0.05,
//            100000,
//            1000) << std::endl;
//    Experiments::PredPreyAssimilation();
    //    Experiments::GnuplotTest();
//    Experiments::PredPreyExpt();
//    Experiments::CatMouseExpt();
//    Experiments::RandomWalk();

    return 0;
}

