//
// Created by daniel on 15/03/2022.
//

#ifndef ABMCMC_TRAJECTORYIMPORTANCE_H
#define ABMCMC_TRAJECTORYIMPORTANCE_H

#include "Trajectory.h"
#include "IndexSet.h"
#include "ABM.h"

template<class AGENT>
class TrajectoryImportance {

    Trajectory<AGENT>            eventTrajectory;
    StateTrajectory<AGENT>       stateTrajectory;
    std::vector<double>                     stateLogProbs;     // factors in the log probability by agent time/state
    double                                  totalLogProb;

    IndexSet                                statesToUpdate;
    IndexSet                                updatedActs;

    TrajectoryImportance(const std::vector<ABM::occupation_type> &eventTrajectory):
        eventTrajectory(eventTrajectory),
        stateTrajectory(eventTrajectory) { }

    // log of P(X)/P_i(X)
    double logImportance() {

    }

    // register a perturbation to the state
    void perturb(const SparseVec<ABM::occupation_type> &changes) {

    }

};


#endif //ABMCMC_TRAJECTORYIMPORTANCE_H
