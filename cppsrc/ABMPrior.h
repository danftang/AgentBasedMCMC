//
// Created by daniel on 27/07/2021.
//

#ifndef GLPKTEST_ABMPRIOR_H
#define GLPKTEST_ABMPRIOR_H

#include <math.h>
#include "ConvexPMF.h"
#include "Trajectory.h"
#include "ABMConstraints.h"

// Prior PMF over ABM trajectories, given a PMF over start states
template<typename AGENT>
class ABMPrior {
public:

    static ConvexPMF PMF(ConvexPMF startStatePMF, int nTimesteps) {
        return ConvexPMF([startState = std::move(startStatePMF.logProb)](const std::vector<double> &X) {
            const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
            return T.logProb() + startState(T(0));
            },
                         nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()+1,
                         ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps) +
                                 startStateConstraintsToTrajectoryConstraints(startStatePMF.convexSupport));
    }

    static std::function<std::vector<double>()> sampler(std::function<std::vector<double>()> startStateSampler, int nTimesteps) {
        return [startSampler = std::move(startStateSampler),nTimesteps]() {
            return Trajectory<AGENT>(nTimesteps, ModelState<AGENT>(startSampler()));
        };
    }


    static glp::Constraint startStateConstraintToTrajectoryConstraint(const glp::Constraint &startStateConstraint) {
        glp::LinearSum trajectoryCoeffs;
        for(int i=0; i<startStateConstraint.coefficients.sparseSize(); ++i) {
            trajectoryCoeffs += startStateConstraint.coefficients.values[i]*State<AGENT>(0, startStateConstraint.coefficients.indices[i]);
        }
        return startStateConstraint.lowerBound <= trajectoryCoeffs <= startStateConstraint.upperBound;
    }

    static ConvexPolyhedron startStateConstraintsToTrajectoryConstraints(const ConvexPolyhedron &startStateConstraints) {
        ConvexPolyhedron trajectoryConstraints;
        for(const glp::Constraint &constraint: startStateConstraints) {
            trajectoryConstraints.push_back(startStateConstraintToTrajectoryConstraint(constraint));
        }
        return trajectoryConstraints;
    }
};


#endif //GLPKTEST_ABMPRIOR_H
