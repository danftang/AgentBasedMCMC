//
// Created by daniel on 27/07/2021.
//

#ifndef GLPKTEST_ABMPRIOR_H
#define GLPKTEST_ABMPRIOR_H

#include <math.h>
#include "ConvexPMF.h"
#include "Trajectory.h"

// Prior PMF over ABM trajectories, given a PMF over start states
template<typename AGENT, typename STARTPMF>
class ABMPrior: public ConvexPMF {
public:
    STARTPMF startStatePMF;
    int nTimesteps;

    ABMPrior(STARTPMF startState, int nTimesteps):
    ConvexPMF([&](const std::vector<double> &X) { return this->trajectoryLogProb(reinterpret_cast<const Trajectory<AGENT> &>(X)); }),
    startStatePMF(std::move(startState)),
    nTimesteps(nTimesteps) {
        convexSupport.addConstraints(continuityConstraints(nTimesteps));
        convexSupport.addConstraints(interactionConstraints(nTimesteps));
        convexSupport.addConstraints(actFermionicConstraints(nTimesteps));
    }

    std::vector<double> nextSample() {
        return Trajectory<AGENT>(nTimesteps, startStatePMF.nextSample());
    }


    double trajectoryLogProb(const Trajectory<AGENT> &T) {
        assert(T.nTimesteps() == nTimesteps);
        const double infeasibilityPenalty = 1.0;
        StateTrajectory<AGENT> stateTrajectory(T);
        double logP = startStatePMF.logProb(stateTrajectory[0]);
        for(int t=0; t<nTimesteps; ++t) {
            for(int agentId = 0; agentId<AGENT::domainSize(); ++agentId) {
                std::vector<double> actPMF = AGENT(agentId).timestep(stateTrajectory[t], infeasibilityPenalty);
                for(int actId=0; actId<AGENT::actDomainSize(); ++actId) {
                    double occupation = fabs(T[Event<AGENT>(t, agentId, actId)]);
                    if(occupation > tol) {
                        logP += occupation * log(actPMF[actId]);
//                   - CombinatoricsUtils.factorialLog(X[eventId]) // add this for non-act-fermionic trajectories
                    }
                }
            }
        }

        // If any state occupation number, m, is greater than 1 then we need to
        // multiply the prob by !m since the m agents can be assigned to m acts in
        // !m ways.
        for(const ModelState<AGENT> &step: stateTrajectory) {
            for(double occupation: step) {
                if(fabs(occupation) > 1.0 + tol) {
                    logP += lgamma(fabs(occupation) + 1.0);
                }
            }
        }

        return logP;
    }

protected:

    static std::vector<glp::Constraint> continuityConstraints(int nTimesteps) {
        std::vector<glp::Constraint> constraints;
        std::vector<std::vector<int>> incomingEdges = consequencesByEndState();
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                glp::Constraint &constraint = constraints.emplace_back(0.0 ,0.0);
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    constraint.coefficients.insert(Event<AGENT>(time, agentState, act), 1.0);
                }
                // incoming edges
                int timeOffset = (time-1)*AGENT::domainSize()*AGENT::actDomainSize();
                for (int inEdge: incomingEdges[agentState]) {
                    constraint.coefficients.insert(timeOffset + inEdge, -1.0);
                }
            }
        }
        return constraints;
    }

    static std::vector<glp::Constraint> interactionConstraints(int nTimesteps) {
        std::vector<glp::Constraint> constraints;
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    for(const glp::Constraint &actConstraint : agent.constraints(time, act)) {
                        push_xImpliesY(constraints, Event(time,agent,act), actConstraint);
                    }
                }
            }
        }
        return constraints;
    }

    // Returns constraint x -> y
    // under the assumption that
    // 0 <= x <= 1
    // and 0 <= y_i <= 1
    // by using the identity
    //
    static void push_xImpliesY(std::vector<glp::Constraint> &constraints, int x, const glp::Constraint &y) {
        if(y.upperBound != INFINITY) {
            glp::Constraint &upperBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] > 0.0) upperBoundConstraint.upperBound += y.coefficients.values[i];
                upperBoundConstraint.coefficients.insert(y.coefficients.indices[i], y.coefficients.values[i]);
            }
            upperBoundConstraint.coefficients.insert(x, upperBoundConstraint.upperBound - y.upperBound);
        }
        if(y.lowerBound != -INFINITY) {
            glp::Constraint &lowerBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] < 0.0) lowerBoundConstraint.upperBound -= y.coefficients.values[i];
                lowerBoundConstraint.coefficients.insert(y.coefficients.indices[i], -y.coefficients.values[i]);
            }
            lowerBoundConstraint.coefficients.insert(x, lowerBoundConstraint.upperBound + y.lowerBound);
        }
    }


    static std::vector<glp::Constraint> actFermionicConstraints(int nTimesteps) {
        std::vector<glp::Constraint> constraints;
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    Event<AGENT> event = Event<AGENT>(time,agentState,act);
                    constraints.push_back(0.0 <= 1.0*event <= 1.0);
                }
            }
        }
        return constraints;
    }


    static std::vector<std::vector<int>> consequencesByEndState() {
        std::vector<std::vector<int> > endStateToEvents(AGENT::domainSize());
        std::vector<AGENT> consequences;
        for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
            AGENT agent(agentState);
            for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                consequences = agent.consequences(act);
                for(const AGENT &endState: consequences) {
                    endStateToEvents[endState].push_back(Event(0,agent,act));
                }
            }
        }
        return endStateToEvents;
    }


};


#endif //GLPKTEST_ABMPRIOR_H
