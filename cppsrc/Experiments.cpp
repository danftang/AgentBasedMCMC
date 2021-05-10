//
// Created by daniel on 06/05/2021.
//

#include <vector>
#include "Experiments.h"
#include "Observation.h"
#include "CatMouseAgent.h"
#include "ABMProblem.h"
#include "SimplexMCMC.h"

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
    abm.simplex();
    std::cout << "LP relaxation initial trajectory: " << abm.primalSolution() << std::endl;
    abm.intOpt();
    std::cout << "MIP initial trajectory: " << abm.mipSolution() << std::endl;


//    Trajectory<CatMouseAgent> initialTrajectory;
//    initialTrajectory.add(Event(0,leftCat, CatMouseAgent::STAYPUT),1.0);
//    std::cout << "Initial trajectory is" << std::endl;
//    std::cout << initialTrajectory << std::endl;

//    abm.stdBasis(); // TODO: make this automatic
//    abm.warmUp();

//    SimplexMCMC mcmc(abm, initialTrajectory);
//
//    for(n in 1..6) {
//        val sample = mcmc.nextSample()
//        CatAndMouseABM.plot(sample)
//        println(sample)
//    }

}
