//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_EVENT_H
#define GLPKTEST_EVENT_H

#include <utility>

template<typename AGENT>
class Event {
    int id;
    Event(int time, const AGENT &agent, const typename AGENT::Act &act):
    id(time*AGENT::domainSize*AGENT::Act::domainSize + agent*AGENT::domainSize + act) { }

    operator int() const { return id; }
};

//template<typename AGENT>
//Event<AGENT> event(int time, const AGENT &agent, const typename AGENT::Act &act) {
//    return Event(time, agent, act);
//}

template<typename AGENT>
std::pair<const double, const Event<AGENT>> operator *(const double &multiplier, const Event<AGENT> &event) {
    return std::pair(multiplier,event);
}

#endif //GLPKTEST_EVENT_H
