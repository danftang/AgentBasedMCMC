//
// Created by daniel on 13/05/2021.
//

#ifndef GLPKTEST_STATETRAJECTORY_H
#define GLPKTEST_STATETRAJECTORY_H

#include <vector>
#include <map>

#include "glpkpp.h"
#include "State.h"
#include "ModelState.h"
#include "SimplexMCMC.h"

// A vector of timesteps, whith each timestep a map from agent state to occupation number
template<typename AGENT>
class StateTrajectory: public std::vector<ModelState<AGENT>> {
public:
    using std::vector<ModelState<AGENT>>::operator [];

//    static constexpr double tol = SimplexMCMC::tol;

//    StateTrajectory(const glp::SparseVec &actTrajectory) {
//        for(int i=1; i <= actTrajectory.sparseSize(); ++i) {
//            auto event = Event<AGENT>(actTrajectory.indices[i]);
//            if(event.time() >= this->size()) this->resize(event.time()+1);
//            (*this)[event.time()][event.agent()] += actTrajectory.values[i];
//        }
//    }


    StateTrajectory(const std::vector<double> &eventTrajectory):
    std::vector<ModelState<AGENT>>((eventTrajectory.size() - 1) / (AGENT::domainSize() * AGENT::actDomainSize())) {
        for(int eventId=1; eventId < eventTrajectory.size(); ++eventId) {
            double occupation = fabs(eventTrajectory[eventId]);
            if(occupation > tol) {
                auto event = Event<AGENT>(eventId);
                (*this)[event.time()][event.agent()] += occupation;
            }
        }
    }


    double operator [](const State<AGENT> &agentState) const {
        if(agentState.time >= this->size()) return 0.0;
        return (*this)[agentState.time][agentState.agent];
    }

    double &operator [](const State<AGENT> &agentState) {
        return (*this)[agentState.time][agentState.agent];
    }
};


#endif //GLPKTEST_STATETRAJECTORY_H
