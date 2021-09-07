//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_EVENT_H
#define GLPKTEST_EVENT_H

#include <utility>

template<typename AGENT>
class Event: public glp::X {
public:
    Event(int time, const AGENT agent, const typename AGENT::Act act):
    X(time*AGENT::domainSize()*AGENT::actDomainSize() + agent*AGENT::actDomainSize() + (int)act + 1) { }

    Event(int eventId): X(eventId) {}

    int time() const        { return (id-1)/(AGENT::domainSize()*AGENT::actDomainSize()); }
    AGENT agent() const     { return ((id-1)%(AGENT::domainSize()*AGENT::actDomainSize()))/AGENT::actDomainSize(); }
    typename AGENT::Act act() const { return (id-1)%AGENT::actDomainSize(); }
    std::vector<AGENT> consequences() const { return agent().consequences(act()); }

    friend std::ostream &operator <<(std::ostream &out, const Event &event) {
        out << event.agent() << ":" << event.act() << "@T" << event.time();
        return out;
    }
};


#endif //GLPKTEST_EVENT_H
