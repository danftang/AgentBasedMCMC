// An example of an agent definition
//
// A Domain is a type that
// - has a non-explicit integer constructor
// - is implicitly convertible to int.

// An agent should be a Domain that includes
// AGENT::Act - type of act that this agent performs, should be a Domain
//  void consequences(AGENT::Act act, std::vector<AGENT> consequences); // consequences of an action
// int domainSize - number of Agent states
// int actDomainSize - number of act states
//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_CATMOUSEAGENT_H
#define GLPKTEST_CATMOUSEAGENT_H


#include <vector>
#include <set>
#include "../ModelState.h"
#include "../version1/Constraint.h"
#include "../Trajectory.h"

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
    static constexpr int domainSize= 4;
    static constexpr int actDomainSize = 2;

    CatMouseAgent(int ordinal): stateId(ordinal) {}
    CatMouseAgent(Type type, Position position): stateId(position + 2*type) { }
    operator int() const { return stateId; }
    Position position() const { return Position(stateId%2); }
    Type type() const { return Type(stateId/2); }

//    std::vector<double> timestep(const ModelState<CatMouseAgent> &others) const;
//    std::vector<double> timestep(const Trajectory<CatMouseAgent> &others, int time) const;
//    double marginalTimestep(Act act) const;
    std::vector<CatMouseAgent> consequences(Act act) const; // the consequences of an act
//    // returns the constraints implied by the given act
//    std::vector<Constraint<ABM::occupation_type>> constraints(int time, Act act) const; // to be generated automatically by static analysis...eventually.
//
//   std::vector<CatMouseAgent> neighbours() const;

    template<class TRAJECTORY>
    static double logEventProb(const Event<CatMouseAgent> &event, const TRAJECTORY &trajectory);

    template<class DOMAIN>
    static std::vector<int> eventProbDependencies(const Event<CatMouseAgent> &event);

    friend std::ostream &operator <<(std::ostream &out, const CatMouseAgent &agent) {
        out << (agent.type()==CAT?"CAT":"MSE") << ":" << (agent.position()==LEFT?"L":"R");
        return out;
    }

};


template<class TRAJECTORY>
double CatMouseAgent::logEventProb(const Event<CatMouseAgent> &event, const TRAJECTORY &trajectory) {
    if (event.agent().type() == CAT) return (event.act() == MOVE)?log(pCatMove):log(1.0-pCatMove);
    if (trajectory[State(event.time(), CatMouseAgent(CAT, event.agent().position()))] >= 1.0) {
        return (event.act() == MOVE)?0.0:-INFINITY;
    }
    return (event.act() == MOVE)?-INFINITY:0.0;
}

template<class DOMAIN>
std::vector<int> CatMouseAgent::eventProbDependencies(const Event<CatMouseAgent> &event) {
    if (event.agent().type() == CAT) return {};
    return {DOMAIN::indexOf(State(event.time(), CatMouseAgent(CAT, event.agent().position())))};
}


#endif //GLPKTEST_CATMOUSEAGENT_H
