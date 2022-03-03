//
// Created by daniel on 28/07/2021.
//

#ifndef GLPKTEST_BINOMIALAGENT_H
#define GLPKTEST_BINOMIALAGENT_H


#include <vector>
#include "glpkpp.h"
#include "../ModelState.h"

// Agent that has no interactions and can either stay put or move right (on a circular, 1D domain)
class BinomialAgent {
public:
    static int GRIDSIZE;
    static double pMove; // probability of moving right

    typedef int Act;

    // Agent Domain stuff
    static int domainSize() { return GRIDSIZE; }
    static constexpr int actDomainSize() { return 2; }

    int stateId;

    BinomialAgent(int id): stateId(id) {}

    operator int() const { return stateId; }

    std::vector<double> timestep(const ModelState<BinomialAgent> &others) const {
        return {1.0-pMove, pMove};
    }

    std::vector<double> marginalTimestep() const {
        return {1.0-pMove, pMove};
    }

    std::vector<BinomialAgent> consequences(Act act) const {
        if(act == 1) return {BinomialAgent((stateId + 1)%GRIDSIZE)};
        return { *this };
    }

    std::vector<glp::Constraint> constraints(int time, Act act) const {
        return { };
    }

//    friend std::ostream &operator <<(std::ostream &out, const BinomialAgent &binomialAgent) {
//        out << binomialAgent.stateId;
//    }
};


#endif //GLPKTEST_BINOMIALAGENT_H
