//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_OBSERVATION_H
#define GLPKTEST_OBSERVATION_H

#include <boost/math/distributions/binomial.hpp>
#include "SimplexMCMC.h"
#include "Event.h"

template<typename AGENT>
class Observation {
public:
    int time;
    AGENT agent;
    int numberObserved;
    double pObserve;

    Observation(int time, const AGENT &agent, int numberObserved, const double &pObserve):
    time(time),
    agent(agent),
    numberObserved(numberObserved),
    pObserve(pObserve) {}

    Observation(int time, const AGENT &agent, const Trajectory<AGENT> &trajectory) {

    }

//    val footprintsObserved: Boolean = (Random.nextDouble() < pObserveFootprints(time, lookedFor, trajectory))
//
//    override fun logLikelihood(trajectory: Trajectory<PredPreyAgent, Acts>): Double {
//        return if(footprintsObserved) {
//            ln(pObserveFootprints(time, lookedFor, trajectory))
//        } else {
//            ln(1.0-pObserveFootprints(time, lookedFor, trajectory))
//        }
//    }
//
//    override fun eventConstraints(): List<MutableConstraint<Fraction>> {
//        return if(footprintsObserved) {
//            listOf(MutableConstraint(hashMapOf(lookedFor.ordinal to Fraction.ONE), ">=", Fraction.ONE).stateToEventConstraint(time))
//        } else {
//            emptyList()
//        }
//    }
//
//    companion object {
//            val pObserveIfPresent = 0.9
//
//            fun pObserveFootprints(time: Int, lookedFor: PredPreyAgent, trajectory: Trajectory<PredPreyAgent, Acts>): Double {
//                return if(trajectory.nAgents(time, lookedFor) >= 1) pObserveIfPresent else 0.0
//            }
//    }
//

    double logLikelihood(const glp::SparseVec &trajectory) const {
        int n = trajectory[State(time,agent)];
        if(n < numberObserved) return -std::numeric_limits<double>::infinity();
        return log(boost::math::pdf(boost::math::binomial(n, pObserve), numberObserved));
    }
    
    std::vector<glp::Constraint> constraints() const;
};


#endif //GLPKTEST_OBSERVATION_H
