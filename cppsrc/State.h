//
// Created by daniel on 04/05/2021.
//

#ifndef GLPKTEST_STATE_H
#define GLPKTEST_STATE_H

#include <utility>
#include "LinearSum.h"
#include "Event.h"

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
        int beginIndex = Event<AGENT>(time, agent, 0);
        for(int actId=0; actId < AGENT::actDomainSize(); ++actId) {
            sum.add(beginIndex+actId, c);
        }
        return sum;
    }

    friend glp::LinearSum operator *(double c, const State<AGENT> &state) {
        return state * c;
    }

    friend std::ostream &operator <<(std::ostream &out, const State &state) {
        out << state.agent << "@t" << state.time;
        return out;
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
