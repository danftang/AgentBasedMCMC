//
// Created by daniel on 29/03/2022.
//

#ifndef ABMCMC_AGENTDATAFLOW_H
#define ABMCMC_AGENTDATAFLOW_H

#include "Dataflow.h"
#include "../ABM.h"
#include "../ModelState.h"

namespace dataflow {

    template<class AGENT>
    class TrajectoryToModelState: public Transform {
    public:
        int time;
        int nTimesteps;

        TrajectoryToModelState(int NTimesteps, int Time): time(Time), nTimesteps(NTimesteps) {
            assert(Time <= NTimesteps && Time >= 0);
        }

        // consumer doesn't have to take r-value, compiler will do whatever conversion is necessary.
        template<class CONSUMER>
        bool operator()(CONSUMER &modelStateConsumer, const std::vector<ABM::occupation_type> &trajectory) {
            return modelStateConsumer(ModelState<AGENT>(trajectory, nTimesteps, time));
        }

    };
}

#endif //ABMCMC_AGENTDATAFLOW_H
