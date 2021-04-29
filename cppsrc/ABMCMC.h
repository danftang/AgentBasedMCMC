//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKTEST_ABMCMC_H
#define GLPKTEST_ABMCMC_H


#include "Observation.h"
#include "SimplexMCMC.h"

template<typename AGENT>
class ABMCMC {

    std::vector<Observation<AGENT>> observations;
    glpkpp::GlpProblem  problem;
    SimplexMCMC         simplex;

    int nSamples = 0;
    int nRejections = 0;
    int fractionalRunLength = 0;
    int nTimesteps;

    ABMCMC(int nTimesteps, std::vector<Observation<AGENT>> &observations,glpkpp::SparseVec &initialSample) :
    problem(ABMProblem(nTimesteps,observations)),
    simplex(problem),
    nTimesteps(nTimesteps) {

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
    static glpkpp::GlpProblem ABMProblem(int nTimesteps, const std::vector<Observation<AGENT>> &observations) {
        glpkpp::GlpProblem prob;
        prob.ensureNVars(AGENT::domainSize * AGENT::Act::domainSize * nTimesteps);
        addContinuityConstraints(prob, nTimesteps);
        addInteractionConstraints(prob, nTimesteps);
        return prob;
    }

    static void addContinuityConstraints(glpkpp::GlpProblem &problem, int nTimesteps) {
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

    static void addInteractionConstraints(glpkpp::GlpProblem &problem, int nTimesteps) {
        for(int time = 0; time < nTimesteps; ++time) {
            for (int agentState = 0; agentState < AGENT::domainSize; ++agentState) {
                AGENT agent(agentState);
                for (int act = 0; act < AGENT::Act::domainSize; ++act) {
                    for(const glpkpp::Constraint &actConstraint : agent.constraints(act)) {
                        addXImpliesY(problem, eventId(time,agentState,act), actConstraint);
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
    static void addXImpliesY(glpkpp::GlpProblem &problem, int x, const glpkpp::Constraint &y) {
            glpkpp::Constraint upperBoundConstraint(-std::numeric_limits<double>::infinity(), 0.0);
            glpkpp::Constraint lowerBoundConstraint(0.0, std::numeric_limits<double>::infinity());
            for(auto term: y.coefficients) {
                if (term.second > 0.0)
                    upperBoundConstraint.upperBound += term.second;
                else
                    lowerBoundConstraint.lowerBound += term.second;
                upperBoundConstraint.coefficients[term.first] = term.second;
                lowerBoundConstraint.coefficients[term.first] = -term.second;
            }
            upperBoundConstraint.coefficients[x] = upperBoundConstraint.upperBound - y.upperBound;
            lowerBoundConstraint.coefficients[x] = lowerBoundConstraint.lowerBound + y.lowerBound;
            problem.addConstraint(upperBoundConstraint);
            problem.addConstraint(lowerBoundConstraint);
    }



    static void calcConsequencesByEndState(std::vector<std::vector<int> > &endStateToEvents) {
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


    static int eventId(int time, int state, int act) {
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
