//
// Created by daniel on 27/07/2021.
//

#ifndef GLPKTEST_ABMWINDOW_H
#define GLPKTEST_ABMWINDOW_H

#include <math.h>
#include "ConvexPMF.h"
#include "Trajectory.h"
#include "ABMConstraints.h"
#include "AgentStateObservation.h"

// Represents the (measurable) space of trajectories of an ABM with a given number of timesteps,
// from which we can generate various useful PMFs
template<typename AGENT>
class ABMWindow {
public:
    int nTimesteps;

    ABMWindow(int nTimesteps): nTimesteps(nTimesteps) {
        assert(nTimesteps != 0); // don't allow null space
    }


    ConvexPMF prior(ConvexPMF startStatePMF) {
        return ConvexPMF([startState = std::move(startStatePMF.logProb)](const std::vector<double> &X) {
            const Trajectory<AGENT> &T = reinterpret_cast<const Trajectory<AGENT> &>(X);
            return T.logProb() + startState(T(0));
            },
                         dimension(),
                         ABMConstraints<AGENT>::actFermionicABMConstraints(nTimesteps) +
                         ABMConstraints<AGENT>::startStateConstraintsToTrajectoryConstraints(startStatePMF.convexSupport));
    }


    std::function<std::vector<double>()> priorSampler(std::function<std::vector<double>()> startStateSampler) {
        return [startSampler = std::move(startStateSampler),nTimesteps = this->nTimesteps]() {
            return Trajectory<AGENT>(nTimesteps, ModelState<AGENT>(startSampler()));
        };
    }


    ConvexPMF likelihood(const AgentStateObservation<AGENT> &observation) {
        ConvexPMF::PMF logP;
        if(observation.state.time == nTimesteps) {
            logP = [observation](const std::vector<double> &X) {
                return observation.logP(observation.state.backwardOccupationNumber(X));
            };
        } else {
            logP = [observation](const std::vector<double> &X) {
                return observation.logP(observation.state.forwardOccupationNumber(X));
            };
        }
        return ConvexPMF(std::move(logP), dimension(), observation.support());
    }


    int dimension() {
        return nTimesteps*AGENT::domainSize()*AGENT::actDomainSize()+1;
    }


};


#endif //GLPKTEST_ABMWINDOW_H
