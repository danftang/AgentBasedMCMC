//
// Created by daniel on 18/05/2021.
//

#ifndef GLPKTEST_PREDPREYAGENT_H
#define GLPKTEST_PREDPREYAGENT_H

#include <vector>
#include <map>
#include "glpkpp.h"
#include "../ModelState.h"

class PredPreyAgent {
public:
    static int GRIDSIZE;

    typedef int Act;


    enum ActNames {
        DIE,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        GIVEBIRTH
    };


    enum Type {
        PREDATOR,
        PREY
    };

    int stateId;

    // Agent Domain stuff
    static int domainSize() { return GRIDSIZE*GRIDSIZE*2; }
    static constexpr int actDomainSize() { return 6; }

    PredPreyAgent(int ordinal): stateId(ordinal) {}
    PredPreyAgent( int x, int y, Type type): stateId(x + GRIDSIZE*y + GRIDSIZE*GRIDSIZE*type + 1) { }
    operator int() const { return stateId; }
    int xPosition() const { return (stateId-1)%GRIDSIZE; }
    int yPosition() const { return ((stateId-1)/GRIDSIZE)%GRIDSIZE; }
    Type type() const { return Type((stateId-1)/(GRIDSIZE*GRIDSIZE)); }

    int xLeft() const { return(((stateId-1)+GRIDSIZE-1)%GRIDSIZE); }
    int xRight() const { return(stateId%GRIDSIZE); }
    int yUp() const { return(((stateId-1)/GRIDSIZE+1)%GRIDSIZE); }
    int yDown() const { return(((stateId-1)/GRIDSIZE + GRIDSIZE-1)%GRIDSIZE); }


    std::vector<double> timestep(const ModelState<PredPreyAgent> &others) const;
    std::vector<PredPreyAgent> consequences(Act act) const; // the consequences of an act
    // returns the constraints implied by the given act
    std::vector<glp::Constraint> constraints(int time, Act act) const; // to be generated automatically by static analysis...eventually.

    friend std::ostream &operator <<(std::ostream &out, const PredPreyAgent &agent) {
        out << agent.type() << ":" << agent.xPosition() << "," << agent.yPosition();
        return out;
    }


    double surroundingCountOf(Type type, const ModelState<PredPreyAgent> &others) const {
        return others[PredPreyAgent(xRight(),yPosition(),type)] +
        others[PredPreyAgent(xLeft(),yPosition(),type)] +
        others[PredPreyAgent(xPosition(),yUp(),type)] +
        others[PredPreyAgent(xPosition(),yDown(),type)];
    }

//
//
//    // returns a constraint in terms of state occupation numbers
//    // Manually generated for now...


};


#endif //GLPKTEST_PREDPREYAGENT_H
