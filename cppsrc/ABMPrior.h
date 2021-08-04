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
template<typename AGENT, typename STARTSTATEPMF>
class ABMPrior: public ConvexPMF {
public:
    typedef ABMPrior<AGENT,STARTSTATEPMF> DefaultSampler;

    STARTSTATEPMF startStatePMF; // a sampleable, convex PMF
    int nTimesteps;


    ABMPrior(STARTSTATEPMF startStateDistribution, int nTimesteps):
    ConvexPMF([this](const std::vector<double> &X) { return (*this)(X); },
              ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps)),
    startStatePMF(std::move(startStateDistribution)),
    nTimesteps(nTimesteps) {
        for(const glp::Constraint &constraint: startStatePMF.convexSupport) {
            convexSupport.push_back(startStateConstraintToTrajectoryConstraint(constraint));
        }
    }


    ABMPrior(const ModelState<AGENT> &fixedStartState, int nTimesteps):
            ConvexPMF(this->trajectoryLogProb),
            startStatePMF(),
            nTimesteps(nTimesteps) {
    }



    std::vector<double> nextSample() {
        return Trajectory<AGENT>(nTimesteps, ModelState<AGENT>(startStatePMF.nextSample()));
    }


    double operator()(const std::vector<double> &X) {
        const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
        assert(T.nTimesteps() == nTimesteps);
        return T.logProb() + startStatePMF(T(0));
    }


    static glp::Constraint startStateConstraintToTrajectoryConstraint(const glp::Constraint &startStateConstraint) {
        glp::LinearSum trajectoryCoeffs;
        for(int i=0; i<startStateConstraint.coefficients.sparseSize(); ++i) {
            trajectoryCoeffs += startStateConstraint.coefficients.values[i]*State<AGENT>(0, startStateConstraint.coefficients.indices[i]);
        }
        return startStateConstraint.lowerBound <= trajectoryCoeffs <= startStateConstraint.upperBound;
    }

};


#endif //GLPKTEST_ABMPRIOR_H
