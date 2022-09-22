//
// Created by daniel on 28/07/2021.
//

#ifndef GLPKTEST_BINOMIALAGENT_H
#define GLPKTEST_BINOMIALAGENT_H


#include <vector>
#include "../ModelState.h"
#include "../Constraint.h"
#include "../ABM.h"
#include "../TrajectoryDependencies.h"

// Agent that has no interactions and can either stay put or move right (on a circular, 1D domain)
template<int GRIDSIZE>
class BinomialAgent {
public:
    static const double pMove; // probability of moving right
    static const double lpMove; // log probability of moving right
    static const double lpNotMove; // log probability of not moving


    typedef int Act;

    // Agent Domain stuff
    static int domainSize() { return GRIDSIZE; }
    static constexpr int actDomainSize() { return 2; }

    int stateId;

    BinomialAgent(int id): stateId(id) {}

    operator int() const { return stateId; }

//    std::vector<double> timestep(const ModelState<BinomialAgent> &others) const {
//        return {1.0-pMove, pMove};
//    }

    template<class TRAJECTORY>
    static double logEventProb(const Event<BinomialAgent<GRIDSIZE>> &event, const TRAJECTORY &trajectory) {
        return event.act() == 1?lpMove:lpNotMove;
    }

    static TrajectoryDependencies<BinomialAgent<GRIDSIZE>> eventProbDependencies(const Event<BinomialAgent<GRIDSIZE>> &event) {
        return {{},{}};
    }

    std::vector<BinomialAgent<GRIDSIZE>> consequences(Act act) const {
        if(act == 1) return { BinomialAgent<GRIDSIZE>((stateId + 1)%GRIDSIZE) };
        return { *this };
    }

};

template<int GRIDSIZE>
const double BinomialAgent<GRIDSIZE>::pMove = 0.5;
template<int GRIDSIZE>
const double BinomialAgent<GRIDSIZE>::lpMove = log(pMove);
template<int GRIDSIZE>
const double BinomialAgent<GRIDSIZE>::lpNotMove = log(1.0-pMove);

#endif //GLPKTEST_BINOMIALAGENT_H
