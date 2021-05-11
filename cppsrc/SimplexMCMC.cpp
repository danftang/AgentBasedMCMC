//
// Created by daniel on 15/04/2021.
//

#include "SimplexMCMC.h"
#include "PivotPoint.h"
#include "Random.h"
#include <algorithm>
#include <cmath>

SimplexMCMC::SimplexMCMC(glp::Problem &prob) :
        Simplex(prob),
        colLastNonZero(n,-1),
        rowLatestCompletionPivot(m, -1),
        latestCompletionBegin(m, -1),
        lnRowPivotCount(m,0.0),
        probability(*this) {
    for(int i = 1; i<=m; ++i) {
        currentBasis[head[i]] = i;
    }
}


void SimplexMCMC::randomPivot() {

}

// Calculate the degeneracy probability of the current state
// from a cold start (i.e. nothing already calculated)
double SimplexMCMC::lnDegeneracyProb() {
//    double row[n+1];
    int nextVarBegin;
    int nextVarEnd = n+m+1;
    int i,j,k,last,lastj;
    int c = 0;
    double lnP = 0.0;

    auto basisVar = currentBasis.rbegin();
    while(basisVar != currentBasis.rend()) {
        i = basisVar->second;
        std::vector<double> row = tableauRow(i);
        ++basisVar;
        nextVarBegin = basisVar->first + 1;
        last = 0;
        for(j = 1; j <=n; ++j) {
            if(row[j] != 0.0) {
                k = head[m+j];
                if(colLastNonZero[j] == -1) {
                    colLastNonZero[j] = i;
                    if(k >= nextVarBegin && k < nextVarEnd) ++c;
                }
                if(k>last) {
                    last = k;
                    lastj = j;
                }
            }
        }
        rowLatestCompletionPivot[i] = last;
        lnP += lnRowPivotCount[i] = std::log(c);
        if(last < nextVarEnd) nextVarEnd = last;
        latestCompletionBegin[i] = nextVarEnd;
        pivot(i,lastj,false);
    }
    return lnP;
}


// Choose a pivot
// and reject based on Metropolis-Hastings
// returns the next sample
glp::SparseVec SimplexMCMC::nextSample() {
    PivotPoint proposalPivot = proposePivot();
    processProposal(proposalPivot);
    ++nSamples;
    if(probability.logFractionalPenalty != 0.0) {
        ++fractionalRunLength;
        std::cout << "${nSamples} In fractional state. Fractional penalty ${state.logFractionalPenalty}";
    }
    return probability.currentSample;
}


void SimplexMCMC::processProposal(PivotPoint proposalPivot) {
    double acceptanceDenominator = probability.logProb() + log(transitionProb(proposalPivot));
    double revTransitionProb = reverseTransitionProb(proposalPivot);
    pivot(proposalPivot.i, proposalPivot.j, false);
    BasisProbability destinationProb(*this);
    double acceptanceNumerator = destinationProb.logProb() + log(revTransitionProb);
    double logAcceptance = std::min(acceptanceNumerator - acceptanceDenominator, 0.0);
//    if(isnan( logAcceptance )) println("NaN Acceptance $acceptanceNumerator / $acceptanceDenominator logPiv = ${state.logProbOfPivotState} transition prob = ${transitionProb(revertState.reversePivot)} columnWeight = ${columnWeights.P(revertState.reversePivot.col)} nPivots = ${nPivots(revertState.reversePivot.col)}");
    if (!(logAcceptance == NAN) && Random::nextDouble() >= exp(logAcceptance)) { // explicity accept if both numerator and denominator are -infinity
        // reject
//        if(revertState.cache.logFractionalPenalty != 0.0 && state.logFractionalPenalty == 0.0) {
//            println("Rejecting fraction->integer transition with denominator = ${revertState.cache.logPX} + ${revertState.cache.logDegeneracyProb} + ${revertState.cache.logFractionalPenalty} + ${acceptanceDenominator - revertState.cache.logFractionalPenalty - revertState.cache.logDegeneracyProb - revertState.cache.logPX} = $acceptanceDenominator, numerator = ${state.logPX} + ${state.logDegeneracyProb} + ${ln(transitionProb(revertState.reversePivot))} = $acceptanceNumerator, acceptance = $logAcceptance")
//        }
//        if(revertState.cache.logFractionalPenalty == 0.0 && state.logFractionalPenalty != 0.0) {
//            println("Rejecting integer->fraction transition with denominator = ${revertState.cache.logPX} + ${revertState.cache.logDegeneracyProb} + ${acceptanceDenominator - revertState.cache.logFractionalPenalty - revertState.cache.logDegeneracyProb - revertState.cache.logPX} = $acceptanceDenominator, numerator = ${state.logPX} + ${state.logDegeneracyProb} +  + ${revertState.cache.logFractionalPenalty} + ${ln(transitionProb(revertState.reversePivot))} = $acceptanceNumerator, acceptance = $logAcceptance")
//        }
//        revertLastPivot()
    } else {
        // Accept proposal
    }
}

double SimplexMCMC::transitionProb(PivotPoint proposal) {
    return 0;
}

double SimplexMCMC::reverseTransitionProb(PivotPoint proposal) {
    return 0;
}

PivotPoint SimplexMCMC::proposePivot() {
    // choose a pivot column
    PivotPoint pivot(0, Random::nextInt(nVars() - nRows()));
    std::vector<double> pivotCol = tableauCol(pivot.j);

    return PivotPoint(0, 0);
}


std::vector<int> SimplexMCMC::calcPivotRows(int j, const std::vector<double> &colVec) {
    std::vector<int> pivRows;
    int k = head[m+j];
    double db_dDSign = isAtLowerBound(j)?-1.0:1.0;
    double deltaMin = u[k] - l[k];
    for(int i=1; i<=nRows();++i) {
        double db_dD = colVec[i] * db_dDSign;
        if(db_dD != 0.0) {
            double Db = db_dD < 0.0 ? (l[head[i]] - b[i]) : (u[head[i]] - b[i]);
            double dD = Db / db_dD;
            if (dD < deltaMin) {
                pivRows.clear();
                pivRows.push_back(i);
                deltaMin = dD;
            } else if (dD == deltaMin) {
                pivRows.push_back(i);
            }
        }
    }
    return pivRows;
}