// Represents the set of binary-valued solutions to a
// set of constraints.
// DOMAIN should be a type that has a subscript operator
// that can be assigned integers
// SUPPORT should be a type with isValidSolution(DOMAIN)
//
// Created by daniel on 10/08/2021.
//

#ifndef GLPKTEST_BINARYSOLUTIONSET_H
#define GLPKTEST_BINARYSOLUTIONSET_H

#include <bitset>
#include <cassert>
#include "EqualityConstraints.h"

template<class DOMAIN, class SUPPORT>
class BinarySolutionSet {
public:
    class Iterator {
    public:
        unsigned long solutionId;
        DOMAIN solution;
        const SUPPORT &support;


        Iterator(const DOMAIN &zeroState, const SUPPORT &support, int id):
        solutionId(id),
        solution(zeroState),
        support(support)
        {
            assert(solution.size() < 32);
            recalculateTrajectory();
//            std::cout << "Made iterator with support:\n" << support << std::endl;
        }

        // construct as begin()
        Iterator(const DOMAIN &zeroState, const SUPPORT &support):
        Iterator(zeroState, support, -1)
        {
            // move forward to first valid exactEndState (or end)
            ++(*this);
        }

        const DOMAIN &operator *() { return solution; }

        Iterator &operator ++() {
            do {
                ++solutionId;
                recalculateTrajectory();
            } while(solutionId < (1 << solution.size()) && !support.isValidSolution(solution));
            return *this;
        }



        bool operator !=(const Iterator &other) const {
            return solutionId != other.solutionId;
        }


        void recalculateTrajectory() {
            std::bitset<32> trajBits(solutionId);
            for(int d=0; d < solution.size(); ++d) solution[d] = trajBits[d];
        }
    };


    const SUPPORT &support;
    DOMAIN zeroState;

    BinarySolutionSet(const SUPPORT &Support, DOMAIN zeroState): support(Support), zeroState(zeroState) {
    }

    Iterator begin() { return Iterator(zeroState, support); }
    Iterator end() { return Iterator(zeroState, support, 1 << zeroState.size());  }

};


#endif //GLPKTEST_BINARYSOLUTIONSET_H
