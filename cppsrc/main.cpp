#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"
#include "Experiments.h"
#include "StlStream.h"
#include "Random.h"

using glp::X;

int main() {
//    Experiments::CatMouseExpt();
    Experiments::Pivot();

    return 0;
}

