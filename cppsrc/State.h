//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_STATE_H
#define GLPKTEST_STATE_H

#include <utility>
#include "LinearSum.h"

template<typename AGENT>
class State {
public:
    int       time;
    AGENT     agent;

    State(int time, const AGENT &agent):
    time(time),
    agent(agent) { }

    glp::LinearSum operator *(double c) const {
        glp::LinearSum sum;
        int beginIndex = Event(time, agent, typename AGENT::Act(0));
        for(int actId=0; actId < (int)AGENT::Act::domainSize; ++actId) {
            sum.add(beginIndex+actId, c);
        }
        return sum;
    }

    friend glp::LinearSum operator *(double c, const State<AGENT> &state) {
        return state * c;
    }


//    friend LinearSum &operator +=(LinearSum &sum, const std::pair<double,State<AGENT>> &stateTerm) {
//        int beginIndex = Event(stateTerm.second.time, stateTerm.second.agent, 0);
//        for(int actId=0; actId < AGENT::Act::domainSize; ++actId) {
//            sum += LinearSum(stateTerm.first, beginIndex+actId);
//        }
//        return sum;
//    }

};



#endif //GLPKTEST_STATE_H
