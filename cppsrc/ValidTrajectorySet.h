//
// Created by daniel on 10/08/2021.
//

#ifndef GLPKTEST_VALIDTRAJECTORYSET_H
#define GLPKTEST_VALIDTRAJECTORYSET_H

#include <bitset>
#include "AssimilationProblem.h"

class ValidTrajectorySet {
public:
    class TrajectoryIterator {
    public:
        unsigned long trajectoryId;
        std::vector<double> trajectory;
        const ConvexPolyhedron &support;


        TrajectoryIterator(int nDimensions, const ConvexPolyhedron &support, int id)
        : trajectoryId(id),
        trajectory(nDimensions),
        support(std::move(support)) {
            assert(nDimensions <= 32);
            trajectory[0] = 0.0;
            recalculateTrajectory();
            std::cout << "Made iterator with support:\n" << support << std::endl;
        }

        TrajectoryIterator(int nDimensions, const ConvexPolyhedron &support)
        : TrajectoryIterator(nDimensions, support, 0) {
            assert(nDimensions <= 32);
        }

        const std::vector<double> &operator *() { return trajectory; }

        TrajectoryIterator &operator ++() {
            do {
                ++trajectoryId;
                recalculateTrajectory();
//                std::cout << "Trying trajectory " << trajectory << std::endl;
            } while(trajectoryId < (1<<(trajectory.size()-1)) && !support.isValidSolution(trajectory));
            return *this;
        }


//        bool operator ==(const TrajectoryIterator &other) const {
//            return trajectoryId == other.trajectoryId;
//        }

        bool operator !=(const TrajectoryIterator &other) const {
            return trajectoryId != other.trajectoryId;
        }


        void recalculateTrajectory() {
            std::bitset<32> trajBits(trajectoryId);
            for(int d=1; d<trajectory.size(); ++d) trajectory[d] = trajBits[d-1];
        }
    };


    const AssimilationProblem &problem;
    ConvexPolyhedron support;

    ValidTrajectorySet(const AssimilationProblem &problem)
    : problem(problem),
    support(problem.PMF().convexSupport) {
    }

    TrajectoryIterator begin() { return TrajectoryIterator(problem.dimension(), support); }
    TrajectoryIterator end() {return TrajectoryIterator(problem.dimension(), support, 1<<problem.dimension());  }

};


#endif //GLPKTEST_VALIDTRAJECTORYSET_H
