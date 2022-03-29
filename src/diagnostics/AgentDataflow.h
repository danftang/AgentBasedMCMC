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
        TrajectoryToModelState(int Time): time(Time) { }

        // consumer doesn't have to take r-value, compiler will do whatever conversion is necessary.
        bool operator()(std::function<bool(ModelState<AGENT> &&)> modelStateConsumer, const std::vector<ABM::occupation_type> &trajectory) {
            return modelStateConsumer(ModelState<AGENT>(trajectory, time));
        }

    };
}

#endif //ABMCMC_AGENTDATAFLOW_H
