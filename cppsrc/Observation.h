//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_OBSERVATION_H
#define GLPKTEST_OBSERVATION_H

#include <boost/math/distributions/binomial.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "SimplexMCMC.h"
#include "Event.h"
#include "Trajectory.h"
#include "State.h"

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

    // make observation of the trajectory
    Observation(int time, const AGENT &agent,const double &pObserve, const Trajectory<AGENT> &trajectory):
    time(time),
    agent(agent),
    pObserve(pObserve)
    {
        boost::random::binomial_distribution binom(trajectory(time,agent),pObserve);
        numberObserved = binom(boost::mt11213b());
    }

//
//    override fun eventConstraints(): List<MutableConstraint<Fraction>> {
//        return if(footprintsObserved) {
//            listOf(MutableConstraint(hashMapOf(lookedFor.ordinal to Fraction.ONE), ">=", Fraction.ONE).stateToEventConstraint(time))
//        } else {
//            emptyList()
//        }
//    }
//

    double logLikelihood(const glp::SparseVec &trajectory) const {
        int n = trajectory[State(time,agent)];
        if(n < numberObserved) return -std::numeric_limits<double>::infinity();
        return log(boost::math::pdf(boost::math::binomial(n, pObserve), numberObserved));
    }

    // if we observed m agents, there cannot be less than m agents present (though there may be more if
    // pObserve < 1.0)
    std::vector<glp::Constraint> constraints() const {
        return std::vector<glp::Constraint>({ glp::Constraint(
                numberObserved,
                1.0 * State(time, agent),
                std::numeric_limits<double>::infinity()
                ) });
    }
};


#endif //GLPKTEST_OBSERVATION_H
