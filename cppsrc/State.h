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
    static const std::vector<std::vector<Event<AGENT>>> incomingEventsByState;

    int       time;
    AGENT     agent;

    State(int time, const AGENT &agent):
    time(time),
    agent(agent) { }


    double forwardOccupationNumber(const std::vector<double> &trajectory) const;

    double backwardOccupationNumber(const std::vector<double> &trajectory) const;

    double occupationUpperBound() const {
        return std::min(AGENT::actDomainSize()*1.0, incomingEventsByState[agent].size()*1.0);
    }


    glp::LinearSum operator *(double c) const {
        glp::LinearSum eventVector;
        int beginIndex = Event<AGENT>(time, agent, 0);
        for(int actId=0; actId < AGENT::actDomainSize(); ++actId) {
            eventVector += c*glp::X(beginIndex+actId);
        }
        return eventVector;
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


    static std::vector<std::vector<Event<AGENT>>> calculateIncomingEventsByState() {
        std::vector<std::vector<Event<AGENT>>> endStateToEvents(AGENT::domainSize());
        std::vector<AGENT> consequences;
        for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
            AGENT agent(agentState);
            for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                consequences = agent.consequences(act);
                for(const AGENT &endState: consequences) {
                    endStateToEvents[endState].push_back(Event(0,agent,act));
                }
            }
        }
        return endStateToEvents;
    }

};


template<typename AGENT>
const std::vector<std::vector<Event<AGENT>>> State<AGENT>::incomingEventsByState = State<AGENT>::calculateIncomingEventsByState();

template<typename AGENT>
double State<AGENT>::backwardOccupationNumber(const std::vector<double> &trajectory) const {
    assert(time != 0);
    double occupation = 0;
    for(const Event<AGENT> &incomingEvent : incomingEventsByState[agent]) {
        occupation += trajectory[Event(time-1, incomingEvent.agent(), incomingEvent.act())];
    }
    return occupation;
}

template<typename AGENT>
double State<AGENT>::forwardOccupationNumber(const std::vector<double> &trajectory) const {
    double occupation = 0;
    for(int act=0; act<AGENT::actDomainSize(); ++act) {
        occupation += trajectory[Event(time,agent,act)];
    }
    return occupation;
}


#endif //GLPKTEST_STATE_H
