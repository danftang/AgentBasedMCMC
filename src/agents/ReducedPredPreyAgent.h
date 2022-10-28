//
// Created by daniel on 02/10/22.
//

#ifndef ABMCMC_REDUCEDPREDPREYAGENT_H
#define ABMCMC_REDUCEDPREDPREYAGENT_H

#include <vector>
#include <ostream>
#include <cassert>
#include <cmath>
#include "../Event.h"

template<int GRIDSIZE>
class ReducedPredPreyAgent {
public:
    typedef int Act;
    static constexpr int gridsize = GRIDSIZE;

    enum ActNames {
        GIVEBIRTH,
        MOVELEFT,
        MOVERIGHT,
        MOVEUP,
        MOVEDOWN,
        DIE
    };
    static constexpr int actDomainSize= 6;

    enum Type {
        PREDATOR,
        PREY
    };
    static constexpr int typeDomainSize = 2;

    static constexpr double lpMove = -1.79175946922805500081; // ln(1/6)
    static constexpr double lpPredBirth = -1.79175946922805500081;
    static constexpr double lpPreyBirth = -1.79175946922805500081;
    static constexpr double lpPredDeath = -1.79175946922805500081;
    static constexpr double lpPreyDeath = -1.79175946922805500081;

    static constexpr int domainSize = GRIDSIZE*GRIDSIZE*2;

    int stateId;

    ReducedPredPreyAgent(int ordinal): stateId(ordinal) {}
    ReducedPredPreyAgent(int x, int y, Type type): stateId(x + GRIDSIZE * y + GRIDSIZE * GRIDSIZE * type) { }

    operator int() const { return stateId; }

    std::vector<ReducedPredPreyAgent<GRIDSIZE>> consequences(Act act) const; // the consequences of an act

    template<class TRAJECTORY>
    static double logEventProb(const Event<ReducedPredPreyAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory);

    template<class DOMAIN>
    static std::vector<int> eventProbDependencies(const Event<ReducedPredPreyAgent<GRIDSIZE>> &event) { return {}; }

    friend std::ostream &operator <<(std::ostream &out, const ReducedPredPreyAgent<GRIDSIZE> &agent) {
        out << (agent.type()==PREY?"PREY":"PRED") << ":(" << agent.xPosition() << "," << agent.yPosition() <<")";
        return out;
    }

    int xPosition() const { return stateId%GRIDSIZE; }
    int yPosition() const { return (stateId/GRIDSIZE)%GRIDSIZE; }
    Type type() const { return Type(stateId/(GRIDSIZE*GRIDSIZE)); }
    Type otherType() const { return type()==PREDATOR?PREY:PREDATOR; }

    int xLeft() const { return((stateId+GRIDSIZE-1)%GRIDSIZE); }
    int xRight() const { return((stateId+1)%GRIDSIZE); }
    int yUp() const { return((stateId/GRIDSIZE+1)%GRIDSIZE); }
    int yDown() const { return((stateId/GRIDSIZE + GRIDSIZE-1)%GRIDSIZE); }
};


template<int GRIDSIZE>
template<class TRAJECTORY>
double ReducedPredPreyAgent<GRIDSIZE>::logEventProb(const Event<ReducedPredPreyAgent<GRIDSIZE>> &event,
                                                    const TRAJECTORY &trajectory) {
    switch (event.act()) {
        case MOVEUP:
        case MOVEDOWN:
        case MOVELEFT:
        case MOVERIGHT:
            return lpMove;
        case GIVEBIRTH:
            return event.agent().type() == PREDATOR? lpPredBirth : lpPreyBirth;
        case DIE:
            return event.agent().type() == PREDATOR ? lpPredDeath : lpPreyDeath;
    }
    assert(false);
    return -INFINITY;
}


template<int GRIDSIZE>
std::vector<ReducedPredPreyAgent<GRIDSIZE>> ReducedPredPreyAgent<GRIDSIZE>::consequences(Act act) const {
    switch (act) {
        case DIE:
            return std::vector<ReducedPredPreyAgent<GRIDSIZE>>();
        case MOVELEFT:
            return std::vector<ReducedPredPreyAgent<GRIDSIZE>>(
                    {ReducedPredPreyAgent<GRIDSIZE>(xLeft(), yPosition(), type())});
        case MOVERIGHT:
            return std::vector<ReducedPredPreyAgent<GRIDSIZE>>(
                    {ReducedPredPreyAgent<GRIDSIZE>(xRight(), yPosition(), type())});
        case MOVEUP:
            return std::vector<ReducedPredPreyAgent<GRIDSIZE>>(
                    {ReducedPredPreyAgent<GRIDSIZE>(xPosition(), yUp(), type())});
        case MOVEDOWN:
            return std::vector<ReducedPredPreyAgent<GRIDSIZE>>(
                    {ReducedPredPreyAgent<GRIDSIZE>(xPosition(), yDown(), type())});
        case GIVEBIRTH:
            return std::vector<ReducedPredPreyAgent<GRIDSIZE>>(
                    {*this, ReducedPredPreyAgent<GRIDSIZE>(xRight(), yPosition(), type())});
    }
    assert(false); // unrecognized act
    return std::vector<ReducedPredPreyAgent<GRIDSIZE>>();
}


#endif //ABMCMC_REDUCEDPREDPREYAGENT_H
