//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_EVENT_H
#define GLPKTEST_EVENT_H

#include <utility>

template<typename AGENT>
class Event: public glp::X {
public:
    Event(int time, const AGENT &agent, const typename AGENT::Act &act):
    X(time*AGENT::domainSize*(int)AGENT::Act::domainSize + agent*AGENT::domainSize + (int)act) { }


};

#endif //GLPKTEST_EVENT_H
