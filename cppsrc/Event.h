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
    X(time*AGENT::domainSize()*AGENT::actDomainSize() + agent*AGENT::domainSize() + (int)act + 1) { }

    Event(int eventId): X(eventId) {}

    int time() const { return id/(AGENT::domainSize()*AGENT::actDomainSize()); }
    AGENT agent() const { return (id%(AGENT::domainSize()*AGENT::actDomainSize()))/AGENT::actDomainSize(); }
    typename AGENT::Act act() const { return id%AGENT::actDomainSize(); }

    friend std::ostream &operator <<(std::ostream &out, const Event &event) {
        out << "[" << event.time() << ", " << event.agent() << " -> " << event.act() << "]";
        return out;
    }
};


#endif //GLPKTEST_EVENT_H
