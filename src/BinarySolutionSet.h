// Represents the set of binary-valued solutions to a ConvexPolyhedron
//
// Created by daniel on 10/08/2021.
//

#ifndef GLPKTEST_BINARYSOLUTIONSET_H
#define GLPKTEST_BINARYSOLUTIONSET_H

#include <bitset>
#include <cassert>
#include "ConvexPolyhedron.h"

template<class T>
class BinarySolutionSet {
public:
    class Iterator {
    public:
        unsigned long solutionId;
        std::vector<T> solution;
        const ConvexPolyhedron<T> &support;


        Iterator(int nDimensions, const ConvexPolyhedron<T> &support, int id)
        : solutionId(id),
          solution(nDimensions),
          support(std::move(support)) {
            assert(nDimensions < 32);
            recalculateTrajectory();
//            std::cout << "Made iterator with support:\n" << support << std::endl;
        }

        // construct as begin()
        Iterator(int nDimensions, const ConvexPolyhedron<T> &support)
        : Iterator(nDimensions, support, -1) {
            // move forward to first valid exactEndState (or end)
            ++(*this);
        }

        const std::vector<T> &operator *() { return solution; }

        Iterator &operator ++() {
            do {
                ++solutionId;
                recalculateTrajectory();
//                std::cout << "Trying solution " << solution << std::endl;
            } while(solutionId < (1 << solution.size()) && !support.isValidSolution(solution));
            return *this;
        }


//        bool operator ==(const Iterator &other) const {
//            return solutionId == other.solutionId;
//        }

        bool operator !=(const Iterator &other) const {
            return solutionId != other.solutionId;
        }


        void recalculateTrajectory() {
            std::bitset<32> trajBits(solutionId);
            for(int d=0; d < solution.size(); ++d) solution[d] = trajBits[d];
        }
    };


    const ConvexPolyhedron<T> &support;
    int nDimensions;

    BinarySolutionSet(const ConvexPolyhedron<T> &Support, int dimension = -1): support(Support) {
        if(dimension > 0) nDimensions = dimension; else nDimensions = support.dimension();
    }

    Iterator begin() { return Iterator(nDimensions, support); }
    Iterator end() { return Iterator(nDimensions, support, 1 << nDimensions);  }

};


#endif //GLPKTEST_BINARYSOLUTIONSET_H
