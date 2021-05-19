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
#include "StateTrajectory.h"

template<typename AGENT>
class Observation {
public:
    int time;
    AGENT agent;
    int numberObserved;
    double pObserveIfPresent;

    Observation(int time, const AGENT &agent, int numberObserved, const double &pObserveIfPresent):
            time(time),
            agent(agent),
            numberObserved(numberObserved),
            pObserveIfPresent(pObserveIfPresent) {}

    // make observation of the trajectory
    Observation(int time, const AGENT &agent, double pObserveIfPresent, const Trajectory<AGENT> &trajectory):
            time(time),
            agent(agent),
            pObserveIfPresent(pObserveIfPresent)
    {
        numberObserved = Random::nextBinomial(trajectory(time,agent), pObserveIfPresent);
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

    double logLikelihood(const StateTrajectory<AGENT> &trajectory) const {
        int n = trajectory[time][agent];
        if(n < numberObserved) return -std::numeric_limits<double>::infinity();
        return log(boost::math::pdf(boost::math::binomial(n, pObserveIfPresent), numberObserved));
    }

    // if we observed m agents, there cannot be less than m agents present (though there may be more if
    // pObserve < 1.0)
    std::vector<glp::Constraint> constraints() const {
        return std::vector({ 1.0*State(time,agent) >= numberObserved });
    }

    static std::pair<std::vector<Observation<AGENT>>,Trajectory<AGENT>> generateObservations(const ModelState<AGENT> &startState, int nTimesteps, double pMakeObservation);
};


template<typename AGENT>
std::pair<std::vector<Observation<AGENT>>,Trajectory<AGENT>>
Observation<AGENT>::generateObservations(const ModelState<AGENT> &startState, int nTimesteps, double pMakeObservation) {
        Trajectory<AGENT> trajectory = Trajectory<AGENT>::run(startState, nTimesteps);
        std::vector<Observation<AGENT>> observations;
        observations.reserve(nTimesteps * AGENT::domainSize() * pMakeObservation * 1.5);
        for (int t=0; t<nTimesteps;++t) {
            for (int agentId=1; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    observations.push_back(Observation(t,AGENT(agentId), 0.9, trajectory));
                }
            }
        }
//        checkTrajectorySatisfiesObervations(trajectory, observations)
    return std::make_pair(observations, trajectory);
}


#endif //GLPKTEST_OBSERVATION_H
