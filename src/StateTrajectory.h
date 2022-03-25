//
// Created by daniel on 13/05/2021.
//

#ifndef GLPKTEST_STATETRAJECTORY_H
#define GLPKTEST_STATETRAJECTORY_H

#include <vector>
#include <map>

#include "State.h"
#include "ModelState.h"
#include "ABM.h"

// A vector of timesteps, whith each timestep a map from agent state to occupation number
template<typename AGENT>
class StateTrajectory: public std::vector<ModelState<AGENT>> {
public:
    using std::vector<ModelState<AGENT>>::operator [];


    StateTrajectory(const std::vector<ABM::occupation_type> &eventTrajectory):
    std::vector<ModelState<AGENT>>((eventTrajectory.size() - 1) / (AGENT::domainSize() * AGENT::actDomainSize())) {
        for(int eventId=1; eventId < eventTrajectory.size(); ++eventId) {
            ABM::occupation_type occupation = eventTrajectory[eventId];
            if(occupation != 0) {
                auto event = Event<AGENT>(eventId);
                (*this)[event.time()][event.agent()] += occupation;
            }
        }
    }

    StateTrajectory(int nTimesteps): std::vector<ModelState<AGENT>>(nTimesteps) { }


    ABM::occupation_type operator [](const State<AGENT> &agentState) const {
        if(agentState.time >= this->size()) return 0;
        return (*this)[agentState.time][agentState.agent];
    }

    ABM::occupation_type &operator [](const State<AGENT> &agentState) {
        return (*this)[agentState.time][agentState.agent];
    }

};


#endif //GLPKTEST_STATETRAJECTORY_H
