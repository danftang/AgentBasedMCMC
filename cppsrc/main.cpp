#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"

using glp::X;



int main() {

//    Experiments::GnuplotTest();
    Experiments::PredPreyExpt();
//    Experiments::CatMouseExpt();
//    Experiments::RandomWalk();

    return 0;
}

