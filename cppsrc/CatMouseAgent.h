//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_CATMOUSEAGENT_H
#define GLPKTEST_CATMOUSEAGENT_H


#include <vector>
#include <set>
#include "glpkpp.h"

class CatMouseAgent {
public:
    typedef int Act;

    enum ActNames {
        MOVE,
        STAYPUT
    };

    enum Position {
        LEFT,
        RIGHT
    };

    enum Type {
        CAT,
        MOUSE
    };

    int stateId;

    // Agent Domain stuff
    static constexpr int domainSize() { return 4; }
    static constexpr int actDomainSize() { return 2; }

    CatMouseAgent(int ordinal): stateId(ordinal) {}
    CatMouseAgent(Type type, Position position): stateId(position + 2*type) { }
    operator int() const { return stateId; }
    Position position() const { return Position(stateId%2); }
    Type type() const { return Type(stateId/2); }

    std::vector<double> timestep(std::multiset<CatMouseAgent> others);
    std::vector<CatMouseAgent> consequences(Act act); // the consequences of an act
    // returns the constraints implied by the given act
    std::vector<glp::Constraint> constraints(int time, Act act); // to be generated automatically by static analysis...eventually.

    friend std::ostream &operator <<(std::ostream &out, const CatMouseAgent &agent) {
        out << agent.type() << ":" << agent.position();
        return out;
    }

};



#endif //GLPKTEST_CATMOUSEAGENT_H
