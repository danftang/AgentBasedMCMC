//
// Created by daniel on 30/07/2021.
//

#ifndef GLPKTEST_ABMSTATEDELTAPMF_H
#define GLPKTEST_ABMSTATEDELTAPMF_H

#include "ConvexPMF.h"
#include "ModelState.h"
#include "State.h"

// Represents a PMF that has uniform probability over
// trajectories that have a given model state at a given time, with
// zero probability at all other times.
//
// N.B. The PMF is not normalised.
template<typename AGENT>
class ABMStateDeltaPMF: public ConvexPMF {
public:
    ABMStateDeltaPMF(const ModelState<AGENT> &modelState, int time, int nTimesteps):
    ConvexPMF([](const std::vector<double> &X) { return 0.0; }) {
        convexSupport.ensureNVars(nTimesteps*AGENT::domainSize()*AGENT::actDomainSize());
        for(int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            int occupation = modelState[agentId];
            convexSupport.addConstraint(occupation <= 1.0*State<AGENT>(time, agentId) <= occupation);
        }
    }
};


#endif //GLPKTEST_ABMSTATEDELTAPMF_H
