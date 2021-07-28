//
// Created by daniel on 27/07/2021.
//

#ifndef GLPKTEST_ABMPRIOR_H
#define GLPKTEST_ABMPRIOR_H

#include <math.h>
#include "ConvexPMF.h"
#include "Trajectory.h"

// Prior PMF over ABM trajectories, given a PMF over start states
template<typename AGENT>
class ABMPrior: public ConvexPMF {
public:
//    STARTSTATEPMF startState; // PMF over start states. Shoud implement logPMF and nextSample.
    PMF *t0PMF;
    int nTimesteps;

    template<typename STARTSTATEPMF>
    ABMPrior(STARTSTATEPMF &startState, int nTimesteps): nTimesteps(nTimesteps) {
        t0PMF = new STARTSTATEPMF(startState);
        logPrior = [this](const std::vector<double> &X) { return t0PMF->logProb(X); };
    }

    template<typename STARTSTATEPMF>
    ABMPrior(STARTSTATEPMF &&startState, int nTimesteps): nTimesteps(nTimesteps) {
        t0PMF = new typename std::remove_reference<STARTSTATEPMF>::type(std::move(startState));
    }

    ~ABMPrior() {
        delete t0PMF;
    }

    std::vector<double> nextSample() {
        return Trajectory<AGENT>(nTimesteps, t0PMF->nextSample());
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
                    constraint.coefficients.add(Event<AGENT>(time, agentState, act), 1.0);
                }
                // incoming edges
                int timeOffset = (time-1)*AGENT::domainSize()*AGENT::actDomainSize();
                for (int inEdge: incomingEdges[agentState]) {
                    constraint.coefficients.add(timeOffset + inEdge, -1.0);
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
                upperBoundConstraint.coefficients.add(y.coefficients.indices[i], y.coefficients.values[i]);
            }
            upperBoundConstraint.coefficients.add(x, upperBoundConstraint.upperBound - y.upperBound);
        }
        if(y.lowerBound != -INFINITY) {
            glp::Constraint &lowerBoundConstraint = constraints.emplace_back(-INFINITY, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] < 0.0) lowerBoundConstraint.upperBound -= y.coefficients.values[i];
                lowerBoundConstraint.coefficients.add(y.coefficients.indices[i], -y.coefficients.values[i]);
            }
            lowerBoundConstraint.coefficients.add(x, lowerBoundConstraint.upperBound + y.lowerBound);
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
