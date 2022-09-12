// Generates a map from all non-zero probability solutions of a
// ConstrainedFactoredDistribution on the unit hypercube to the
// exact probability of the solution.
//
// Created by daniel on 19/08/2021.

#ifndef GLPKTEST_EXACTSOLVER_H
#define GLPKTEST_EXACTSOLVER_H

#include <valarray>
#include "BinarySolutionSet.h"
#include "ModelState.h"
#include "StlStream.h"

template<typename DOMAIN>
class ExactSolver {
public:
    typedef std::remove_cv_t<std::remove_reference_t<DOMAIN>> stripped_domain;
    std::map<stripped_domain,double> pmf;

    ExactSolver(const ConstrainedFactorisedDistribution<DOMAIN> &distribution, const DOMAIN &zeroState) {
        std::cout << "In exact solver with distribution\n" << distribution << std::endl;
        double marginalP = 0.0;
        for(const DOMAIN &solution: BinarySolutionSet(distribution.constraints, zeroState)) {
            double jointP = exp(distribution.logPexact(solution));
            if(jointP > 1e-12) {
                assert(distribution.constraints.isValidSolution(solution));
                pmf[solution] = jointP;
                marginalP += jointP;
            }
        }
        for(auto &entry: pmf) { entry.second /= marginalP; }
    }
};


#endif //GLPKTEST_EXACTSOLVER_H
