#include <iostream>
#include "glpkppinclude/glpkpp.h"
#include "CatMouseAgent.h"
#include "ABMProblem.h"
#include "Trajectory.h"

using glp::X;

glp::SparseVec myFunc() {
    return glp::SparseVec(3);
}

int main() {

//    Trajectory<CatMouseAgent> myTraj(3);
//    myTraj = myFunc();
    glp::SparseVec myVec(3);
    Trajectory<CatMouseAgent> myTraj(myVec);

    Observation<CatMouseAgent> myObs(0, CatMouseAgent(0), 1, 0.9);
    auto constraints = myObs.constraints();

    std::copy(constraints.begin(), constraints.end(), std::ostream_iterator<glp::Constraint>(std::cout, ", "));

    std::cout << myObs.logLikelihood(myVec);

//    std::cout << myTraj;

//    ABMProblem<CatMouseAgent> myProb(4, );

//    glp::Problem myProb;
//
//    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
//    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
//    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
//    myProb.addConstraint(0.0 <= 1.0*X(1));
//    myProb.addConstraint(0.0 <= 1.0*X(2));
//    myProb.addConstraint(0.0 <= 1.0*X(3));
//    myProb.setObjective(10.0*X(1) + 6.0*X(2) + 4.0*X(3));
//    myProb.stdBasis();
//    myProb.warmUp();
//
//    glp::Simplex mySimplex(myProb);
//
//    std::cout << myProb;
//    std::cout << mySimplex << std::endl;
//
//    mySimplex.pivot(3,3,true);
//
//    std::cout << mySimplex << std::endl;

    return 0;
}
