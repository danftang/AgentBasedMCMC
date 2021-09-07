//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_CATMOUSEAGENT_H
#define GLPKTEST_CATMOUSEAGENT_H


#include <vector>
#include <set>
#include "glpkpp.h"
#include "../ModelState.h"

class CatMouseAgent {
public:
    static constexpr double pCatMove = 0.25;

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

    std::vector<double> timestep(const ModelState<CatMouseAgent> &others) const;
    std::vector<double> marginalTimestep() const;
    std::vector<CatMouseAgent> consequences(Act act) const; // the consequences of an act
    // returns the constraints implied by the given act
    std::vector<glp::Constraint> constraints(int time, Act act) const; // to be generated automatically by static analysis...eventually.

    friend std::ostream &operator <<(std::ostream &out, const CatMouseAgent &agent) {
        out << (agent.type()==CAT?"CAT":"MSE") << ":" << (agent.position()==LEFT?"L":"R");
        return out;
    }

};



#endif //GLPKTEST_CATMOUSEAGENT_H
