//
// Created by daniel on 19/08/2021.
//

#ifndef GLPKTEST_EXACTSOLVER_H
#define GLPKTEST_EXACTSOLVER_H

#include "BinarySolutionSet.h"

template<typename AGENT>
class ExactSolver {
public:
    ModelState<AGENT> solution;

    ExactSolver(const ConvexPMF &pmf) {
        double marginalP = 0.0;
        for(const std::vector<double> &traj: BinarySolutionSet(pmf.convexSupport, pmf.nDimensions)) {
            double jointP = exp(pmf.logP(traj));
            marginalP += jointP;
            //            std::cout << traj << " " << jointP << std::endl;
            ModelState<AGENT> endState = Trajectory<AGENT>(traj).endState();
            endState *= jointP;
            solution += endState;
        }
        solution *= 1.0/marginalP;
    }

    ExactSolver(const AssimilationWindow<AGENT> &window): ExactSolver(window.posterior) {

    }
};


#endif //GLPKTEST_EXACTSOLVER_H
