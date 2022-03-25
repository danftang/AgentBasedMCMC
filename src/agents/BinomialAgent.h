//
// Created by daniel on 28/07/2021.
//

#ifndef GLPKTEST_BINOMIALAGENT_H
#define GLPKTEST_BINOMIALAGENT_H


#include <vector>
#include "../ModelState.h"
#include "../Constraint.h"
#include "../ABM.h"

// Agent that has no interactions and can either stay put or move right (on a circular, 1D domain)
template<int GRIDSIZE>
class BinomialAgent {
public:
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

    double marginalTimestep(Act act) const {
         if(act == 0) return 1.0-pMove;
         return pMove;
    }

    std::vector<BinomialAgent<GRIDSIZE>> neighbours() {
        return std::vector<BinomialAgent<GRIDSIZE>>();
    }


    std::vector<BinomialAgent<GRIDSIZE>> consequences(Act act) const {
        if(act == 1) return {BinomialAgent<GRIDSIZE>((stateId + 1)%GRIDSIZE)};
        return { *this };
    }

    std::vector<Constraint<ABM::occupation_type>> constraints(int time, Act act) const {
        return { };
    }

//    friend std::ostream &operator <<(std::ostream &out, const BinomialAgent &binomialAgent) {
//        out << binomialAgent.stateId;
//    }
};

template<int GRIDSIZE>
double BinomialAgent<GRIDSIZE>::pMove = 0.5;


#endif //GLPKTEST_BINOMIALAGENT_H
