//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_EVENT_H
#define GLPKTEST_EVENT_H

#include <utility>
#include "X.h"

template<typename AGENT>
class Event: public X {
public:
    Event(int time, const AGENT agent, const typename AGENT::Act act):
            X(time*AGENT::domainSize()*AGENT::actDomainSize() + agent*AGENT::actDomainSize() + (int)act) { }

    Event(int eventId): X(eventId) {}

    int time() const        { return id/(AGENT::domainSize()*AGENT::actDomainSize()); }
    AGENT agent() const     { return id%(AGENT::domainSize()*AGENT::actDomainSize())/AGENT::actDomainSize(); }
    typename AGENT::Act act() const { return id%AGENT::actDomainSize(); }
    std::vector<AGENT> consequences() const { return agent().consequences(act()); }

    friend std::ostream &operator <<(std::ostream &out, const Event &event) {
        out << event.agent() << ":" << event.act() << "@T" << event.time();
        return out;
    }
};

template<class AGENT>
inline LinearSum<ABM::occupation_type> operator *(int coefficient, const Event<AGENT> &variable) {
    return LinearSum<ABM::occupation_type>({ {variable.id, ABM::occupation_type(coefficient)} });
}


#endif //GLPKTEST_EVENT_H
