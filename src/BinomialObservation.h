//
// Created by daniel on 22/07/2021.
//

#ifndef GLPKTEST_BINOMIALOBSERVATION_H
#define GLPKTEST_BINOMIALOBSERVATION_H

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
class BinomialObservation {
public:
    State<AGENT> agentState;
    int numberObserved;
    double pObserveIfPresent;

    BinomialObservation(const State<AGENT> &agentState, int numberObserved, double pObserveIfPresent):
            agentState(agentState),
            numberObserved(numberObserved),
            pObserveIfPresent(pObserveIfPresent) {}

    // make observation of the exactEndState
    BinomialObservation(int time, const AGENT &agent, double pObserveIfPresent, const Trajectory<AGENT> &trajectory):
            agentState(time,agent),
            pObserveIfPresent(pObserveIfPresent)
    {
        numberObserved = Random::nextBinomial(trajectory[State<AGENT>(time,agent)], pObserveIfPresent);
    }

    BinomialObservation(int time, const AGENT &agent, double pObserveIfPresent, int realOccupationNumber):
            agentState(time,agent),
            pObserveIfPresent(pObserveIfPresent)
    {
        numberObserved = Random::nextBinomial(realOccupationNumber, pObserveIfPresent);
    }


    double logLikelihood(const StateTrajectory<AGENT> &trajectory, double infeasibilityLogProb) const {
        int n = trajectory[agentState];
        if(n < numberObserved) {
//            std::cout << "Infeasible observation: observed " << numberObserved << " but only " << n << " present. Logprob = " << infeasibilityLogProb << std::endl;
            return infeasibilityLogProb;
        }
        double logP = log(boost::math::pdf(boost::math::binomial(n, pObserveIfPresent), numberObserved));
//        std::cout << "Observation log likelihood = " << logP << std::endl;
        return logP;
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


    static std::vector<Observation<AGENT>>
    generateObservations(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        int nTimesteps = realTrajectory.nTimesteps();
        std::vector<Observation < AGENT>>
        observations;
        observations.reserve((realTrajectory.size() - 1) * pMakeObservation * 1.5);
        for (int t = 0; t < nTimesteps; ++t) {
            for (int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    AGENT agent(agentId);
                    observations.push_back(Observation(t, agent, pObserveIfPresent, realTrajectory(t, agent)));
                }
            }
        }
//        checkTrajectorySatisfiesObervations(exactEndState, observations)
        debug(std::cout << "Generated observations: " << observations << std::endl);
        return observations;
    }
};


#endif //GLPKTEST_BINOMIALOBSERVATION_H
