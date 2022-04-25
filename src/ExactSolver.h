// Generates a map from all solutions of a WeightedFactoredDistribution to the
// exact probability of the solution.
//
// Created by daniel on 19/08/2021.

#ifndef GLPKTEST_EXACTSOLVER_H
#define GLPKTEST_EXACTSOLVER_H

#include <valarray>
#include "BinarySolutionSet.h"
#include "ModelState.h"
#include "WeightedFactoredConvexDistribution.h"
#include "StlStream.h"

template<typename AGENT>
class ExactSolver {
public:
    std::map<std::vector<ABM::occupation_type>,double> pmf;

    ExactSolver(const WeightedFactoredConvexDistribution<ABM::occupation_type> &distribution) {
        double marginalP = 0.0;
        int nTimesteps = distribution.constraints.dimension() / (AGENT::domainSize() * AGENT::actDomainSize());
        for(const std::vector<ABM::occupation_type> &solution: BinarySolutionSet<ABM::occupation_type>(distribution.constraints)) {
            double jointP = distribution.P(solution);
            pmf[solution] = jointP;
            marginalP += jointP;
            ModelState<AGENT> endState(solution, nTimesteps, nTimesteps);
        }
        for(auto it=pmf.begin(); it != pmf.end(); ++it) {
            it->second /= marginalP;
        }
    }
};


#endif //GLPKTEST_EXACTSOLVER_H
