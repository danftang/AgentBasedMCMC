//
// Created by daniel on 15/04/2021.
//

#include "SimplexMCMC.h"
#include "Pivot.h"
#include "Random.h"
#include "ColumnPivot.h"
#include <algorithm>
#include <cmath>

SimplexMCMC::SimplexMCMC(glp::Problem &prob, std::function<double (const std::vector<double> &)> logProb) :
        Simplex(prob),
        colLastNonZero(n,-1),
        rowLatestCompletionPivot(m, -1),
        latestCompletionBegin(m, -1),
        lnRowPivotCount(m,0.0),
        logProbFunc(logProb),
        probability(*this) {
//    for(int i = 1; i<=m; ++i) {
//        orderedBasis[head[i]] = i;
//    }
}



// Calculate the degeneracy probability of the current state
// from a cold start (i.e. nothing already calculated)
//
// If we split the variables into the ones that are on a bound, d, and those that aren't, x,
// then any basis that includes all x variables will have the same solution.
// However, on an integer solution, all variables are on their bound, so any
// pivot state can represent any state!!!
double SimplexMCMC::lnDegeneracyProb() {
    int nextVarBegin;
    int nextVarEnd = n+m+1;
    int i,j,k,last,lastj;
    int c = 0;
    double lnP = 0.0;


//    auto basisVar = orderedBasis.rbegin();
//    while(basisVar != orderedBasis.rend()) {
//        i = basisVar->second;
//        std::vector<double> row = tableauRow(i);
//        ++basisVar;
//        nextVarBegin = basisVar->first + 1;
//        last = 0;
//        for(j = 1; j <=n-m; ++j) {
//            if(row[j] != 0.0) {
//                k = head[m+j];
//                if(colLastNonZero[j] == -1) {
//                    colLastNonZero[j] = i;
//                    if(k >= nextVarBegin && k < nextVarEnd) ++c;
//                }
//                if(k>last) {
//                    last = k;
//                    lastj = j;
//                }
//            }
//        }
//        rowLatestCompletionPivot[i] = last;
//        lnP += lnRowPivotCount[i] = std::log(c);
//        if(last < nextVarEnd) nextVarEnd = last;
//        latestCompletionBegin[i] = nextVarEnd;
//        glp::Simplex::pivot(i, lastj);
//        lpSolutionIsValid(false);
//    }
    return lnP;
}

double SimplexMCMC::lnFractionalPenalty() {
    // fractions must be in the basis, so check b
    const double tol = 1e-8;
    double penalty = 0.0;
    for(int i=1; i<=m; ++i) {
        if(fabs(round(b[i]) - b[i]) > tol) {
            penalty += fractionalK;
        }
    }
    return penalty;
}


// Choose a pivot
// and reject based on Metropolis-Hastings
// returns the next sample
void SimplexMCMC::nextSample() {
    ColumnPivot proposalPivot = proposePivot();
    processProposal(proposalPivot);
    ++nSamples;
    if(probability.logFractionalPenalty != 0.0) {
        ++fractionalRunLength;
        std::cout << "${nSamples} In fractional state. Fractional penalty ${state.logFractionalPenalty}";
    }
}


void SimplexMCMC::processProposal(const ColumnPivot &proposalPivot) {
    std::cout << "Processing proposal " << proposalPivot.i << ", " << proposalPivot.j << std::endl;
    double acceptanceDenominator = probability.logProb() + log(transitionProb(proposalPivot));
    pivot(proposalPivot);
    ColumnPivot reversePivot = proposalPivot.reverse(*this);
    BasisProbability destinationProb(*this);
    double acceptanceNumerator = destinationProb.logProb() + log(transitionProb(reversePivot));
    double logAcceptance = std::min(acceptanceNumerator - acceptanceDenominator, 0.0);
//    if(isnan( logAcceptance )) println("NaN Acceptance $acceptanceNumerator / $acceptanceDenominator logPiv = ${state.logProbOfPivotState} transition prob = ${transitionProb(revertState.reversePivot)} columnWeight = ${columnWeights.P(revertState.reversePivot.col)} nPivots = ${nPivots(revertState.reversePivot.col)}");
    std::cout << "Log acceptance is " << acceptanceNumerator << " - " << acceptanceDenominator << " = " << logAcceptance << std::endl;
    if (std::isnan(logAcceptance) || Random::nextDouble() <= exp(logAcceptance)) { // explicity accept if both numerator and denominator are -infinity
        // Accept proposal
        std::cout << "Accepting" << std::endl;
        probability = destinationProb;
    } else {
        // reject
        std::cout << "Rejecting" << std::endl;
        pivot(reversePivot);
//        if(revertState.cache.logFractionalPenalty != 0.0 && state.logFractionalPenalty == 0.0) {
//            println("Rejecting fraction->integer transition with denominator = ${revertState.cache.logPX} + ${revertState.cache.logDegeneracyProb} + ${revertState.cache.logFractionalPenalty} + ${acceptanceDenominator - revertState.cache.logFractionalPenalty - revertState.cache.logDegeneracyProb - revertState.cache.logPX} = $acceptanceDenominator, numerator = ${state.logPX} + ${state.logDegeneracyProb} + ${ln(transitionProb(revertState.reversePivot))} = $acceptanceNumerator, acceptance = $logAcceptance")
//        }
//        if(revertState.cache.logFractionalPenalty == 0.0 && state.logFractionalPenalty != 0.0) {
//            println("Rejecting integer->fraction transition with denominator = ${revertState.cache.logPX} + ${revertState.cache.logDegeneracyProb} + ${acceptanceDenominator - revertState.cache.logFractionalPenalty - revertState.cache.logDegeneracyProb - revertState.cache.logPX} = $acceptanceDenominator, numerator = ${state.logPX} + ${state.logDegeneracyProb} +  + ${revertState.cache.logFractionalPenalty} + ${ln(transitionProb(revertState.reversePivot))} = $acceptanceNumerator, acceptance = $logAcceptance")
//        }
//        revertLastPivot()
    }
}

double SimplexMCMC::transitionProb(const ColumnPivot &proposal) {
    return 1.0/(n*proposal.pivotRows.size());
}

//double SimplexMCMC::reverseTransitionProb(Pivot proposal) {
//    return 0;
//}

void SimplexMCMC::randomWalk() {
    pivot(proposePivot());
}

ColumnPivot SimplexMCMC::proposePivot() {
    // choose a pivot column
    int j = Random::nextInt(1,n - m + 1);
    ColumnPivot proposal(*this, j);
    if(proposal.pivotRows.size() >0) {
        proposal.i = proposal.pivotRows[Random::nextInt(proposal.pivotRows.size())];
    } else {
        proposal.i = -1; // incoming var goes to opposite limit.
    }
    return proposal;
}



// assumes no unbounded variables
//std::vector<int> SimplexMCMC::calcPivotRows(int j, const std::vector<double> &colVec) {
//    std::vector<int> pivRows;
//    int kIncoming = head[m + j];
//    double deltaMin = u[kIncoming] - l[kIncoming];
//    double DXi;
//    double DXj;
//    for(int i=1; i<=nRows();++i) {
//        if(colVec[i] != 0.0) {
//            int kOutgoing = head[i];
//            bool outgoingToUpperBound = (colVec[i] > 0.0) ^ isAtUpperBound(j);
//            if(outgoingToUpperBound) {
//                DXi = u[kOutgoing] - b[i];
//            } else {
//                DXi = l[kOutgoing] - b[i];
//            }
//            DXj = fabs(DXi/colVec[i]);
//            if (DXj < deltaMin) {
//                pivRows.clear();
//                pivRows.push_back(i);
//                deltaMin = DXj;
//            } else if (DXj == deltaMin) {
//                pivRows.push_back(i);
//            }
//        }
//    }
//    return pivRows;
//}