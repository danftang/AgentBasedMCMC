//
// Created by daniel on 30/04/2021.
//

#ifndef GLPKTEST_ABMPROBLEM_H
#define GLPKTEST_ABMPROBLEM_H

#include <vector>
#include "glpkpp.h"
#include "Observation.h"
#include "Event.h"

template<typename AGENT>
class ABMProblem: glp::Problem {
    static constexpr double infinity = std::numeric_limits<double>::infinity();

    int nTimesteps;
//    std::vector<Observation<AGENT>> observations;


    ABMProblem(int nTimesteps, const std::vector<Observation<AGENT> > &observations): nTimesteps(nTimesteps) {
        ensureNVars(AGENT::domainSize * AGENT::Act::domainSize * nTimesteps);
        addContinuityConstraints();
        addInteractionConstraints();
    }

    //    fun logProb(X: SparseVector<Fraction>): Double {
//        if(X.isInteger()) {
//            val trajectory = Trajectory(model, X)
//            val prior = trajectory.logPrior()
//            val likelihood = observations.sumByDouble { it.logLikelihood(trajectory) }
////        val penalty = logFractionPenalty(X)
//            val logP = prior + likelihood //+ penalty
////        if(penalty < 0.0) println("Got fractional penalty $penalty")
////        println("prior logprob $prior likelihood logprob $likelihood fraction penalty $penalty = $logP")
//            return logP
//        } else {
//            // TODO: Test for dealing with fractional states
//            return X.nonZeroEntries.values.sumByDouble { it.toDouble() * -2.2 } // constant = average logprob per event
//
//        }
//    }


    void addContinuityConstraints() {
        glp::Constraint constraint(0.0,0.0);
        std::vector<std::vector<int>> incomingEdges;
        calcConsequencesByEndState(incomingEdges);
        for(int time = 0; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                // outgoing edges
                if(time < nTimesteps-1) {
                    for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                        constraint[Event(time, AGENT(agentState), act)] = 1.0;
                    }
                }
                if(time > 0) {
                    // incoming edges
                    int timeOffset = time*AGENT::domainSize*AGENT::Act::domainSize;
                    for (int inEdge: incomingEdges[agentState]) {
                        constraint[timeOffset + inEdge] = -1.0;
                    }
                }
                addConstraint(constraint);
                constraint.coefficients.clear();
            }
        }
    }

    void addInteractionConstraints() {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                    for(const glp::Constraint &actConstraint : agent.constraints(act)) {
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
            for (auto term: y.coefficients) {
                if (term.second > 0.0) upperBoundConstraint.upperBound += term.second;
                upperBoundConstraint.coefficients[term.first] = term.second;
            }
            upperBoundConstraint.coefficients[x] = upperBoundConstraint.upperBound - y.upperBound;
            addConstraint(upperBoundConstraint);
        }
        if(y.lowerBound != -infinity) {
            glp::Constraint lowerBoundConstraint(0.0, infinity);
            for (auto term: y.coefficients) {
                if (term.second < 0.0) lowerBoundConstraint.lowerBound += term.second;
                lowerBoundConstraint.coefficients[term.first] = -term.second;
            }
            lowerBoundConstraint.coefficients[x] = lowerBoundConstraint.lowerBound + y.lowerBound;
            addConstraint(lowerBoundConstraint);
        }
    }




    static void calcConsequencesByEndState(std::vector<std::vector<int> > &endStateToEvents) {
        endStateToEvents.clear();
        endStateToEvents.resize(AGENT::domainSize * AGENT::Act::domainSize);
//        AGENT *pAgent;
        std::vector<AGENT> consequences;
        for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
            AGENT agent(agentState);
            for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                agent.consequences(act, consequences);
                for(AGENT endState: consequences) {
                    endStateToEvents[endState].push_back(Event(0,agent,act));
                }
            }
        }
    }
};


#endif //GLPKTEST_ABMPROBLEM_H
