//
// Created by daniel on 19/08/2021.
//

#ifndef GLPKTEST_EXACTSOLVER_H
#define GLPKTEST_EXACTSOLVER_H

#include "BinarySolutionSet.h"

template<typename AGENT>
class ExactSolver {
public:
    ModelState<AGENT> exactEndState;

    ExactSolver(const ConvexPMF<Trajectory<AGENT>> &pmf) {
        double marginalP = 0.0;
        for(const std::vector<double> &sol: BinarySolutionSet(pmf.convexSupport, pmf.nDimensions)) {
            const Trajectory<AGENT> &traj = reinterpret_cast<const Trajectory<AGENT> &>(sol);
            double logP = pmf.logP(traj);
            double jointP = exp(logP);
            marginalP += jointP;
            assert(!isnan(logP));
//            std::cout << traj << " " << jointP << " " << logP << std::endl;
            ModelState<AGENT> endState = traj.endState();
            endState *= jointP;
            exactEndState += endState;
        }
        exactEndState *= 1.0 / marginalP;
    }

//    ExactSolver(const AssimilationWindow<AGENT> &window): ExactSolver(window.posterior) {
//
//    }
};


#endif //GLPKTEST_EXACTSOLVER_H
