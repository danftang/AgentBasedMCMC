//
// Created by daniel on 18/08/2021.
//

#ifndef GLPKTEST_TRAJECTORYPMF_H
#define GLPKTEST_TRAJECTORYPMF_H

#include "ConvexPMF.h"
#include "Trajectory.h"
#include "ABM.h"

//template<typename AGENT>
//class TrajectoryPMF: public ConvexPMF {
//public:
//    TrajectoryPMF(int nTimesteps, ConvexPMF priorStartState)
//    : ConvexPMF([startState = std::move(priorStartState.extendedLogProb)](const std::vector<double> &X) {
//            const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
//            return T.extendedLogProb() + startState(T(0));
//        },
//                Trajectory<AGENT>::dimension(nTimesteps),
//                ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps) +
//                ABMConstraints<AGENT>::startStateConstraintsToTrajectoryConstraints(priorStartState.convexSupport)) {}
//
//
//};


#endif //GLPKTEST_TRAJECTORYPMF_H
