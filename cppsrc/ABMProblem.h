//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_ABMPROBLEM_H
#define GLPKTEST_ABMPROBLEM_H

#include <vector>
#include "glpkpp.h"
#include "Observation.h"
#include "Event.h"
#include "StlStream.h"

template<typename AGENT>
class ABMProblem: public glp::Problem {
public:
    static constexpr double infinity = std::numeric_limits<double>::infinity();

    int nTimesteps;
    std::vector<Observation<AGENT>> observations;



    ABMProblem(int nTimesteps, std::vector<Observation<AGENT> > observations):
    nTimesteps(nTimesteps),
    observations(observations) {
        ensureNVars(AGENT::domainSize() * AGENT::actDomainSize() * nTimesteps);
        addContinuityConstraints();
        addInteractionConstraints();
        addActFermionicConstraints();
        addObservations(observations);
    }

    double logProb(const glp::SparseVec & X) {
            const auto &trajectory = (const Trajectory<AGENT> &)X;
            double logP = trajectory.logPrior();
            for(auto observation: observations) {
                logP += observation.logLikelihood(trajectory);
            }
            return logP;
    }

protected:

    void addContinuityConstraints() {
        glp::Constraint constraint(0.0,0.0);
        std::vector<std::vector<int>> incomingEdges = consequencesByEndState();
        for(int time = 1; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                // outgoing edges
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    constraint.coefficients.add(Event<AGENT>(time, agentState, act), 1.0);
                }
                // incoming edges
                int timeOffset = (time-1)*AGENT::domainSize()*AGENT::actDomainSize();
                for (int inEdge: incomingEdges[agentState]) {
                    constraint.coefficients.add(timeOffset + inEdge, -1.0);
                }
                addConstraint(constraint);
                constraint.coefficients.clear();
            }
        }
    }

    void addInteractionConstraints() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    for(const glp::Constraint &actConstraint : agent.constraints(time, act)) {
                        addXImpliesY(Event(time,agent,act), actConstraint);
                    }
                }
            }
        }
    }

    // Returns constraint x -> y
    // under the assumption that
    // 0 <= x <= 1
    // and 0 <= y_i <= 1
    // by using the identity
    //
    void addXImpliesY(int x, const glp::Constraint &y) {
        if(y.upperBound != infinity) {
            glp::Constraint upperBoundConstraint(-infinity, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] > 0.0) upperBoundConstraint.upperBound += y.coefficients.values[i];
                upperBoundConstraint.coefficients.add(y.coefficients.indices[i], y.coefficients.values[i]);
            }
            upperBoundConstraint.coefficients.add(x, upperBoundConstraint.upperBound - y.upperBound);
            addConstraint(upperBoundConstraint);
        }
        if(y.lowerBound != -infinity) {
            glp::Constraint lowerBoundConstraint(-infinity, 0.0);
            for (int i=0; i < y.coefficients.sparseSize(); ++i) {
                if (y.coefficients.values[i] < 0.0) lowerBoundConstraint.upperBound -= y.coefficients.values[i];
                lowerBoundConstraint.coefficients.add(y.coefficients.indices[i], -y.coefficients.values[i]);
            }
            lowerBoundConstraint.coefficients.add(x, lowerBoundConstraint.upperBound + y.lowerBound);
            addConstraint(lowerBoundConstraint);
        }
    }


    void addActFermionicConstraints() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
                for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                    Event<AGENT> event = Event<AGENT>(time,agentState,act);
                    addConstraint(0.0 <= 1.0*event <= 1.0);
                    setColKind(event, BINARY);
                }
            }
        }
    }

    void addObservations(const std::vector<Observation<AGENT> > &observations) {
        for(const Observation<AGENT> &observation: observations) {
            for(const glp::Constraint &constraint: observation.constraints()) {
                addConstraint(constraint);
            }
        }
    }


    static std::vector<std::vector<int>> consequencesByEndState() {
        std::vector<std::vector<int> > endStateToEvents;
        endStateToEvents.clear();
        endStateToEvents.resize(AGENT::domainSize());
        std::vector<AGENT> consequences;
        for(int agentState = 0; agentState < AGENT::domainSize(); ++agentState) {
            AGENT agent(agentState);
            for (int act = 0; act < AGENT::actDomainSize(); ++act) {
                consequences = agent.consequences(act);
                for(AGENT endState: consequences) {
                    endStateToEvents[endState].push_back(Event(0,agent,act));
                }
            }
        }
        return endStateToEvents;
    }


};


#endif //GLPKTEST_ABMPROBLEM_H
