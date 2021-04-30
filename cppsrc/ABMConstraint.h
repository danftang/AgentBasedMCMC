//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_ABMCONSTRAINT_H
#define GLPKTEST_ABMCONSTRAINT_H

#include "glpkpp.h"
#include "Event.h"

template<typename AGENT>
class ABMConstraint: glp::Constraint {

    ABMConstraint & operator +=(const std::pair<const double, const Event<AGENT>> &term) {
        Event(1,term.second.agent, term.second.act);
        return *this;
    }
};


#endif //GLPKTEST_ABMCONSTRAINT_H
