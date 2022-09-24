//
// Created by daniel on 06/09/22.
//

#ifndef ABMCMC_NOISELESSAGENTSTATEOBSERVATION_H
#define ABMCMC_NOISELESSAGENTSTATEOBSERVATION_H


#include "ABM.h"
#include "State.h"
#include "ConstrainedFactorisedDistribution.h"
#include "Trajectory.h"
//
//template<typename AGENT>
//class NoiselessAgentStateObservation: public ConstrainedFactorisedDistribution<const Trajectory<AGENT> &,ABM::coefficient_type> {
//public:
//    NoiselessAgentStateObservation(const State<AGENT> &state, ABM::coefficient_type nObserved) {
//        this->addConstraint(state == nObserved);
//    }
//
//    friend std::ostream &operator <<(std::ostream &out, const NoiselessAgentStateObservation<AGENT> & observation) {
//        out << observation.state << " n=" << observation.nObserved;
//        return out;
//    }
//
//};


#endif //ABMCMC_NOISELESSAGENTSTATEOBSERVATION_H
