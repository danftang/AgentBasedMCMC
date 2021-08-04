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
#include "debug.h"

template<typename AGENT>
class Observation {
public:
    State<AGENT> agentState;
    int numberObserved;

    Observation(int time, const AGENT &agent, int numberObserved):
    agentState(time,agent),
    numberObserved(numberObserved) { }

    // make observation of the trajectory
//    Observation(int time, const AGENT &agent, const Trajectory<AGENT> &trajectory):
//    agentState(time,agent),
//    numberObserved(trajectory(time,agent)) { }

    double logLikelihood(const StateTrajectory<AGENT> &trajectory, double infeasibilityLogProb) const {
        return trajectory[agentState]==numberObserved?0.0:infeasibilityLogProb;
    }

    // if we observed m agents, there cannot be less than m agents present (though there may be more if
    // pObserve < 1.0)
    std::vector<glp::Constraint> constraints() const {
        return std::vector({ 1.0*agentState == numberObserved });
    }

    friend std::ostream &operator <<(std::ostream &out, const Observation &obs) {
        out << obs.agentState << " == " << obs.numberObserved;
        return out;
    }


    static std::vector<Observation<AGENT>>
    generateObservations(const Trajectory<AGENT> &realTrajectory, double pMakeObservation) {
        int nTimesteps = realTrajectory.nTimesteps();
        std::vector<Observation<AGENT>> observations;
        observations.reserve((realTrajectory.size()-1)*pMakeObservation*1.5);
        for (int t=0; t<nTimesteps;++t) {
            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    AGENT agent(agentId);
                    observations.push_back(Observation(t,agent, realTrajectory(t,agent)));
                }
            }
        }
//        checkTrajectorySatisfiesObervations(trajectory, observations)
        debug(std::cout << "Generated observations: " << observations << std::endl);
        return observations;
    }


};

//double operator[](const std::vector<double> &X, int time, int agentId ) {
//    return 0.0;
//}

#endif //GLPKTEST_OBSERVATION_H
