//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_CATMOUSEAGENT_H
#define GLPKTEST_CATMOUSEAGENT_H


#include <vector>
#include <set>
#include "glpkpp.h"

class CatMouseAgent {
    enum class Act {
        MOVE,
        STAYPUT,
        domainSize = 2
    };

    enum Position {
        LEFT,
        RIGHT
    };


    Position pos;

    // Agent Domain stuff
    static const int domainSize = 2;
    CatMouseAgent(int ordinal) { if(ordinal%2 == 0) pos = LEFT; else pos = RIGHT;}
    operator int() const { return pos; }


    std::vector<double> timestep(std::multiset<CatMouseAgent> others);
    std::vector<CatMouseAgent> consequences(Act act); // the consequences of an act
    // returns the constraints implied by the given act
    std::vector<glp::Constraint> constraints(Act act); // to be generated automatically by static analysis...eventually.

};


#endif //GLPKTEST_CATMOUSEAGENT_H
