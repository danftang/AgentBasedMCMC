//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include "Experiments.h"
#include "agents/CatMouseAgent.h"
#include "ABMProblem.h"
#include "SimplexMCMC.h"
#include "agents/PredPreyAgent.h"


void Experiments::PredPreyExpt() {
    PredPreyAgent::GRIDSIZE = 16;
    int nTimesteps = 8;
    auto startState = ModelState<PredPreyAgent>::randomPoissonState([](const PredPreyAgent &agent) {
        if(agent.type() == PredPreyAgent::PREDATOR) return 0.02;
        return 0.04;
    });
    std::cout << "Start state: " << startState << std::endl;
    auto [observations, realTrajectory] =
            Observation<PredPreyAgent>::generateObservations(startState, nTimesteps, 0.01);
    ABMProblem<PredPreyAgent> abm(nTimesteps, observations);

    std::cout << "Real trajectory: " << glp::SparseVec(realTrajectory) << std::endl;
    std::cout << "Observations: " << observations << std::endl;
//    std::cout << "Linearised ABM:\n" << abm << std::endl;

    abm.solutionBasis(realTrajectory);
    abm.warmUp();



    SimplexMCMC mcmc(abm, abm.logProbFunc());

    // do we have fixed variables in the basis here?
    for(int i=1; i<=mcmc.m; ++i) {
        if(mcmc.l[mcmc.head[i]] == mcmc.u[mcmc.head[i]]) std::cout << "Fixed var in basis " << i << std::endl;
    }
    for(int j=1; j<=mcmc.n-mcmc.m; ++j) {
        if(mcmc.l[mcmc.head[mcmc.m+j]] == mcmc.u[mcmc.head[mcmc.m+j]]) std::cout << "Fixed null var " << j << std::endl;
    }

    std::cout << glp::SparseVec(mcmc.X()) << std::endl;
    for(int n=0; n<99; ++n) {
        mcmc.nextSample();
        std::cout << glp::SparseVec(mcmc.X()) << std::endl;
//        std::cout << "number of fractional pivots = " << mcmc.countFractionalPivCols() << " / " << mcmc.nNonBasic() << std::endl;
        assert(abm.isValidSolution(mcmc.X()));
    }
}


void Experiments::CatMouseExpt() {
    CatMouseAgent leftCat(CatMouseAgent::Type::CAT, CatMouseAgent::Position::LEFT);
    auto observations = std::vector(
            {Observation(
                    1,
                    leftCat,
                    1,
                    0.95)
            });

    ABMProblem<CatMouseAgent> abm(2, observations);

    std::cout << abm << std::endl;

    // calculate initial trajectory
    abm.cpxBasis();
    abm.simplex();
    std::cout << "LP relaxation initial trajectory: " << abm.primalSolution() << std::endl;
    abm.intOpt();
    std::cout << "MIP initial trajectory: " << abm.mipSolution() << std::endl;
    abm.warmUp();

//    Trajectory<CatMouseAgent> initialTrajectory;
//    initialTrajectory.add(Event(0,leftCat, CatMouseAgent::STAYPUT),1.0);
//    std::cout << "Initial trajectory is" << std::endl;
//    std::cout << initialTrajectory << std::endl;
//    abm.stdBasis(); // TODO: make this automatic
//    abm.warmUp();

    SimplexMCMC mcmc(abm, abm.logProbFunc());
    std::cout << mcmc.X() << std::endl;
    for(int n=0; n<100; ++n) {
        mcmc.nextSample();
//        mcmc.randomWalk();
        std::cout << n << "  " << mcmc.X() << std::endl;
        std::cout << "Valid = " << abm.isValidSolution(mcmc.X()) << std::endl;
        assert(abm.isValidSolution(mcmc.X()));
    }

}


void Experiments::RandomWalk() {
    using glp::X;
    glp::Problem myProb;

    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
    myProb.addConstraint(0.0 <= 1.0*X(1));
    myProb.addConstraint(0.0 <= 1.0*X(2));
    myProb.addConstraint(0.0 <= 1.0*X(3));
    myProb.setObjective(10.0*X(1) + 6.0*X(2) + 4.0*X(3));
    myProb.stdBasis();
    myProb.warmUp();
    std::cout << myProb;

    SimplexMCMC myMCMC(myProb, nullPMF);

    std::cout << myMCMC << std::endl;

    for(int sample=0; sample < 5; ++sample) {
        myMCMC.randomWalk();
        std::cout << myMCMC << std::endl;
        std::cout << "Sample is: " << myMCMC.X() << std::endl;
    }

//    mySimplex.pivot(3,3);
//    std::cout << mySimplex << std::endl;

}