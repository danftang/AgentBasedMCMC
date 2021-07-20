//
// Created by daniel on 20/07/2021.
//

#ifndef GLPKTEST_ASSIMILATIONWINDOW_H
#define GLPKTEST_ASSIMILATIONWINDOW_H

#include <vector>
#include "debug.h"

template<typename AGENT>
class AssimilationWindow {
public:
    Trajectory<AGENT>               realTrajectory;
    std::vector<Observation<AGENT>> observations;
    PoissonState<AGENT>             analysis;

    AssimilationWindow(int nTimesteps,
                     const PoissonState<AGENT> &priorStartState,
                     const ModelState<AGENT> &realStartState,
                     double pMakeObservation,
                     double pObserveIfPresent,
                     int nSamples) :
            realTrajectory(nTimesteps, realStartState),
            observations(Observation<AGENT>::generateObservations(realTrajectory, pMakeObservation, pObserveIfPresent)) {
        ABMProblem<AGENT> abm(nTimesteps, observations, [&](const Trajectory<AGENT> &trajectory) {
            return priorStartState.logProb(trajectory(0));
        });
        SimplexMCMC mcmc(abm, abm.logProbFunc());

        Trajectory<AGENT> firstGuessTrajectory(nTimesteps, priorStartState.sample());
        mcmc.setLPState(firstGuessTrajectory);
        mcmc.findFeasibleStartPoint();
        assert(mcmc.abmSanityChecks());

        for(int n=0; n<nSamples; ++n) {
            mcmc.nextSample();
            debug(if(nSamples<1000 || n%1000 == 1) assert(abm.isValidSolution(mcmc.X())));
            const Trajectory<AGENT> &trajectory = reinterpret_cast<const Trajectory<AGENT> &>(mcmc.X());
            analysis += trajectory.endState();
        }
        debug(std::cout
        << "infeasible/feasible: " << mcmc.infeasibleStatistics.nSamples *100.0/mcmc.feasibleStatistics.nSamples << "%" << std::endl
        << "Feasible sample statistics:" << std::endl << mcmc.feasibleStatistics << std::endl
        << "Infeasible sample statistics:" << std::endl << mcmc.infeasibleStatistics << std::endl;
        );
    }

};


#endif //GLPKTEST_ASSIMILATIONWINDOW_H
