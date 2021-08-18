//
// Created by daniel on 18/08/2021.
//

#ifndef GLPKTEST_TRAJECTORYPRIORPMF_H
#define GLPKTEST_TRAJECTORYPRIORPMF_H

#include "ConvexPMF.h"
#include "Trajectory.h"
#include "ABMConstraints.h"

template<typename AGENT>
class TrajectoryPriorPMF: public ConvexPMF {
public:
    TrajectoryPriorPMF(int nTimesteps, ConvexPMF priorStartState)
    : ConvexPMF([startState = std::move(priorStartState.logProb)](const std::vector<double> &X) {
            const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
            return T.logProb() + startState(T(0));
        },
                Trajectory<AGENT>::dimension(nTimesteps),
                ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps) +
                ABMConstraints<AGENT>::startStateConstraintsToTrajectoryConstraints(priorStartState.convexSupport)) {}


};


#endif //GLPKTEST_TRAJECTORYPRIORPMF_H
