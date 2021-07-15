//
// Created by daniel on 15/07/2021.
//

#ifndef GLPKTEST_ASSIMILATIONWINDOW_H
#define GLPKTEST_ASSIMILATIONWINDOW_H

#include <vector>
#include "ABMProblem.h"
#include "Trajectory.h"
#include "PoissonState.h"

template <typename AGENT>
class AssimilationWindow {
public:

    static PoissonState<AGENT> assimilate(
            int nTimeseps,
            std::vector<Observation<AGENT>> observations,
            const PoissonState<AGENT> &prior,
            const Trajectory<AGENT> &initialTrajectory,
            int nSamples) {

        ABMProblem<AGENT> abm(nTimeseps, observations, [&](const Trajectory<AGENT> &trajectory) {
            return prior.logProb(trajectory(0));
        });
        SimplexMCMC mcmc(abm, abm.logProbFunc());
        mcmc.setLPState(initialTrajectory);

        PoissonState<AGENT> finalState;
        for(int n=0; n<nSamples; ++n) {
            mcmc.nextSample();
            if(nSamples<1000 || n%1000 == 1) {
//                std::cout << "Sample " << n << std::endl;
                assert(abm.isValidSolution(mcmc.X()));
//            std::cout << "Sample " << n << " : " << glp::SparseVec(mcmc.X()) << std::endl;
            }
            const Trajectory<AGENT> &trajectory = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
            finalState += trajectory(abm.nTimesteps);
        }
        std::cout << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl;
        std::cout << "Feasible sample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl;
        std::cout << "Infeasible sample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
        return finalState;
    }



};


#endif //GLPKTEST_ASSIMILATIONWINDOW_H
