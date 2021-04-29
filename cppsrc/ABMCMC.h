//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_ABMCMC_H
#define GLPKTEST_ABMCMC_H


#include "Observation.h"
#include "SimplexMCMC.h"

template<typename AGENT>
class ABMCMC {

    SimplexMCMC         simplex;
    glpkpp::GlpProblem  problem;
    std::vector<Observation<AGENT>> observations;

    int nSamples = 0;
    int nRejections = 0;
    int fractionalRunLength = 0;
    int nTimesteps;

    ABMCMC(int nTimesteps, std::vector<Observation<AGENT>> &observations,
    glpkpp::SparseVec &initialSample) {
        this->nTimesteps = nTimesteps;
        problem.ensureNVars(AGENT::domainSize * AGENT::Act::domainSize * nTimesteps);
//        this.observations = observations
//        val constraints = constraints(model, nTimesteps, observations)
//        this.simplex = IntegerSimplexMCMC(
//                constraints,
//                initialSample
//                ?: run {
//                        ORTools.IntegerSolve(constraints)
//                                .withIndex()
//                                .filter { it.value != 0.0 }
//                        .associate { Pair(it.index, Fraction(it.value)) }
//                        .asVector(FractionOperators)
//                },
//                this::logProb
//        )
//        assert(simplex.isFullyPivoted())
//        assert(simplex.isPrimalFeasible())

//        currentFractionalPenalty = logFractionPenalty(simplex.X(true))
    }


protected:
    void addContinuityConstraints() {
        glpkpp::Constraint constraint(0.0,0.0);
        std::vector<std::vector<int>> incomingEdges;
        calcConsequencesByEndState(incomingEdges);
        for(int time = 0; time < nTimesteps; ++time) {
            for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                // outgoing edges
                if(time < nTimesteps-1) {
                    for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                        constraint[eventId(time, agentState, act)] = 1.0;
                    }
                }
                if(time > 0) {
                    // incoming edges
                    int timeOffset = time*AGENT::domainSize*AGENT::Act::domainSize;
                    for (int inEdge: incomingEdges[agentState]) {
                        constraint[timeOffset + inEdge] = -1.0;
                    }
                }
                problem.addConstraint(constraint);
                constraint.coefficients.clear();
            }
        }
    }

    void addInteractionConstraints() {
        glpkpp::Constraint constraint;
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                    glpkpp::Constraint constraint = agent->constraints(act);

                    // TODO: do act implies constraint and add to LP
                }
            }
        }
    }


    void calcConsequencesByEndState(std::vector<std::vector<int> > &endStateToEvents) {
        endStateToEvents.clear();
        endStateToEvents.resize(AGENT::domainSize * AGENT::Act::domainSize);
        AGENT *pAgent;
        std::vector<AGENT> consequences;
        for(int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
            pAgent = new AGENT(agentState);
            for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                pAgent->consequences(act, consequences);
                for(AGENT endState: consequences) {
                    endStateToEvents[endState].push_back(eventId(0,agentState,act));
                }
            }
            delete pAgent;
        }
    }


    int eventId(int time, int state, int act) {
        return time*(AGENT::domainSize*AGENT::Act::domainSize) +
            state*AGENT::Act::domainSize +
            act;
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
//
//
//
//    fun nextSample(): Trajectory<AGENT,ACT> {
//        return Trajectory(model, simplex.nextIntegerSample())
//    }
//
//
//
//
//    fun<R> expectation(nSamples: Int, initialExpectation: R, expectationAccumulator: (SparseVector<Fraction>, R) -> R): R {
//        var e = initialExpectation
//        var oldSample: SparseVector<Fraction>? = null
//        var rejections = 0
//        var lastTime = Instant.now().toEpochMilli()
//        for(s in 1..nSamples) {
//            val newSample = simplex.nextIntegerSample()
//            e = expectationAccumulator(newSample, e)
//            if(oldSample === newSample) ++rejections
//            oldSample = newSample
//            if(s.rem(100) == 0) {
//                val now = Instant.now().toEpochMilli()
//                println("Got to sample $s in ${(now-lastTime)/1000.0}s largest Numerator,Denominator ${largestNumeratorDenominator()}, Size ${simplex.M.rows.sumBy { it.nonZeroEntries.size }}, Degeneracy ${simplex.degeneracy()} logPiv = ${simplex.state.logProbOfPivotState} logPX = ${simplex.state.logPX}, logPDegeneracy = ${simplex.state.logDegeneracyProb} prob per event = ${simplex.state.logPX/newSample.nonZeroEntries.size}")
//                lastTime = now
////                println(simplex.fractionalLogP - simplex.state.logPX - ln(simplex.transitionProb(simplex.proposePivot())))
//            }
//        }
//        println("Rejection ratio = ${rejections.toDouble()/nSamples}")
//        return e
//    }

};



#endif //GLPKTEST_ABMCMC_H
