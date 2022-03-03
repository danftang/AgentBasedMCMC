//
// Created by daniel on 10/08/2021.
//

#ifndef GLPKTEST_BINARYSOLUTIONSET_H
#define GLPKTEST_BINARYSOLUTIONSET_H

#include <bitset>

class BinarySolutionSet {
public:
    class Iterator {
    public:
        unsigned long solutionId;
        std::vector<double> solution;
        const ConvexPolyhedron &support;


        Iterator(int nDimensions, const ConvexPolyhedron &support, int id)
        : solutionId(id),
          solution(nDimensions),
          support(std::move(support)) {
            assert(nDimensions <= 32);
            solution[0] = 0.0;
            recalculateTrajectory();
//            std::cout << "Made iterator with support:\n" << support << std::endl;
        }

        // construct as begin()
        Iterator(int nDimensions, const ConvexPolyhedron &support)
        : Iterator(nDimensions, support, -1) {
            // move forward to first valid exactEndState (or end)
            ++(*this);
        }

        const std::vector<double> &operator *() { return solution; }

        Iterator &operator ++() {
            do {
                ++solutionId;
                recalculateTrajectory();
//                std::cout << "Trying exactEndState " << exactEndState << std::endl;
            } while(solutionId < (1 << (solution.size() - 1)) && !support.isValidSolution(solution));
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
            for(int d=1; d < solution.size(); ++d) solution[d] = trajBits[d - 1];
        }
    };


//    const AssimilationProblem &problem;
    ConvexPolyhedron support;
    int nDimensions;

    BinarySolutionSet(ConvexPolyhedron support, int dimension = -1)
    : support(support) {
        if(dimension > 0) nDimensions = dimension; else nDimensions = support.dimension();
    }

    Iterator begin() { return Iterator(nDimensions, support); }
    Iterator end() { return Iterator(nDimensions, support, 1 << (nDimensions-1));  }

};


#endif //GLPKTEST_BINARYSOLUTIONSET_H
