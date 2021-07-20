#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"

using glp::X;

int main(int argc, char *argv[]) {

//    std::cout << Experiments::informationIncrease(argc, argv) << std::endl;

    Experiments::PredPreyAssimilation();
    //    Experiments::GnuplotTest();
//    Experiments::PredPreyExpt();
//    Experiments::CatMouseExpt();
//    Experiments::RandomWalk();

    return 0;
}

