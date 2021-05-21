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
    State<AGENT> agentState;
//    int time;
//    AGENT agent;
    int numberObserved;
    double pObserveIfPresent;

    Observation(int time, const AGENT &agent, int numberObserved, const double &pObserveIfPresent):
            agentState(time,agent),
            numberObserved(numberObserved),
            pObserveIfPresent(pObserveIfPresent) {}

    // make observation of the trajectory
    Observation(int time, const AGENT &agent, double pObserveIfPresent, const Trajectory<AGENT> &trajectory):
            agentState(time,agent),
            pObserveIfPresent(pObserveIfPresent)
    {
        numberObserved = Random::nextBinomial(trajectory(time,agent), pObserveIfPresent);
    }


    double logLikelihood(const StateTrajectory<AGENT> &trajectory) const {
        int n = trajectory[agentState];
        if(n < numberObserved) return -std::numeric_limits<double>::infinity();
        return log(boost::math::pdf(boost::math::binomial(n, pObserveIfPresent), numberObserved));
    }

    // if we observed m agents, there cannot be less than m agents present (though there may be more if
    // pObserve < 1.0)
    std::vector<glp::Constraint> constraints() const {
        if(numberObserved == 0) return std::vector<glp::Constraint>();
        return std::vector({ 1.0*agentState >= numberObserved });
    }

    friend std::ostream &operator <<(std::ostream &out, const Observation &obs) {
        out << obs.agentState << " -> " << obs.numberObserved;
        return out;
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
        for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
            if (Random::nextDouble() < pMakeObservation) {
                observations.push_back(Observation(t,AGENT(agentId), 0.9, trajectory));
            }
        }
    }
//        checkTrajectorySatisfiesObervations(trajectory, observations)
    return std::pair(observations, trajectory);
}



#endif //GLPKTEST_OBSERVATION_H
