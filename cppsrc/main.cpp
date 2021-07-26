#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"
#include "UnitTests.h"

using glp::X;


class MyClass {
public:
    int i;

    MyClass(int j): i(j) { std::cout << "Initialising" << std::endl; }
    MyClass(const MyClass &other): i(other.i) { std::cout << "Copying" << std::endl;}
    MyClass(MyClass &&other): i(other.i) { std::cout << "Moving" << std::endl;}

};

int main(int argc, char *argv[]) {

    MyClass i(1234);
    std::function<int()> f = [=]() { return i.i+1; };
//    std::function<int()> g = std::move(f);
    f = [g = std::move(f)]() { return g() + 1; };


    std::cout << f() << std::endl;
//    std::cout << g() << std::endl;

//    UnitTests::testActFermionicDistribution();

//    std::cout << Experiments::informationIncrease(
//            8,
//            2,
//            1,
//            0.16,
//            0.32,
//            0.005,
//            100000,
//            1000) << std::endl;
//    Experiments::PredPreyAssimilation();
    //    Experiments::GnuplotTest();
//    Experiments::PredPreyExpt();
//    Experiments::CatMouseExpt();
//    Experiments::RandomWalk();

    return 0;
}

