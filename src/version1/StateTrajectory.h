//
// Created by daniel on 13/05/2021.
//

#ifndef GLPKTEST_STATETRAJECTORY_H
#define GLPKTEST_STATETRAJECTORY_H

#include <vector>
#include <map>

#include "../State.h"
#include "../ModelState.h"
#include "../ABM.h"

// A vector of timesteps, whith each timestep a map from agent state to occupation number
template<typename AGENT>
class StateTrajectory: public std::vector<ModelState<AGENT>> {
public:
    using std::vector<ModelState<AGENT>>::operator [];


    StateTrajectory(const std::vector<ABM::occupation_type> &eventTrajectory, int nTimesteps):
    std::vector<ModelState<AGENT>>(nTimesteps) {
        for(int eventId=0; eventId < nTimesteps*AGENT::domainSize*AGENT::actDomainSize; ++eventId) {
            ABM::occupation_type occupation = eventTrajectory[eventId];
            if(occupation != 0) {
                auto event = Event<AGENT>(eventId);
                (*this)[event.time()][event.agent()] += occupation;
            }
        }
    }

    StateTrajectory(int nTimesteps): std::vector<ModelState<AGENT>>(nTimesteps) { }


    const ABM::occupation_type &operator [](const State<AGENT> &agentState) const {
//        if(agentState.time >= this->size()) return 0;
        assert(agentState.time < this->size());
        return (*this)[agentState.time][agentState.agent];
    }

    ABM::occupation_type &operator [](const State<AGENT> &agentState) {
        return (*this)[agentState.time][agentState.agent];
    }

//    const ABM::occupation_type &operator [](int index) const {
//        int time = stateTrajectoryIndex / AGENT::domainSize;
//        int state = stateTrajectoryIndex % AGENT::domainSize;
//        return (*this)[time][state];
//    }

};


#endif //GLPKTEST_STATETRAJECTORY_H
