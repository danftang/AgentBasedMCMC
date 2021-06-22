//
// Created by daniel on 15/04/2021.
//

#include "SimplexMCMC.h"
#include "ProposalPivot.h"
#include "Random.h"
#include "ColumnPivot.h"
#include <algorithm>
#include <cmath>

SimplexMCMC::SimplexMCMC(glp::Problem &prob, const std::function<double (const std::vector<double> &)> &logProb) :
        Simplex(prob),
        colLastNonZero(n,-1),
        rowLatestCompletionPivot(m, -1),
        latestCompletionBegin(m, -1),
        lnRowPivotCount(m,0.0),
        logProbFunc(logProb) {
    // TODO: ensure that all auxiliary vars are initially in the basis
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
    do {
        ProbabilisticColumnPivot proposalPivot = proposePivot();
        processProposal(proposalPivot);
    } while(!solutionIsPrimaryFeasible());
    ++nSamples;
//    if(probability.logFractionalPenalty != 0.0) {
//        ++fractionalRunLength;
//        std::cout << "${nSamples} In fractional state. Fractional penalty ${state.logFractionalPenalty}";
//    }
}

// TODO: How to assign probabilities to infeasible states?
// - treat negative occupation numbers as anti-agents
// - have linear approximation of trajectory probability?
// - have non-interacting approximation of probability?
void SimplexMCMC::processProposal(const ProposalPivot &proposalPivot) {
    std::cout << "Processing proposal " << proposalPivot.i << ", " << proposalPivot.j
              << " leavingVarToUpperBound = " << proposalPivot.leavingVarToUpperBound
              << " deltaj = " << proposalPivot.deltaj
              << " incoming var = " << (isAtUpperBound(proposalPivot.j)?u[head[nBasic()+proposalPivot.j]]:l[head[nBasic()+proposalPivot.j]]);
    if(proposalPivot.i > 0) {
        std::cout << " b[i] = " << b[proposalPivot.i] << " leaving var limits " << l[head[proposalPivot.i]]
                  << ":" << u[head[proposalPivot.i]];
    }
    std::cout << std::endl;

    double sourceProb = logProbFunc(X());
    std::cout << "Source LP state is: " << glp::SparseVec(lpSolution) << std::endl;
    updateLPSolution(proposalPivot);
    double destinationProb = logProbFunc(lpSolution);
    std::cout << "Destination LP state is: " << glp::SparseVec(lpSolution) << std::endl;
    revertLPSolution(proposalPivot);
    double logAcceptance = destinationProb - sourceProb + proposalPivot.logAcceptanceContribution;
//    if(isnan( logAcceptance )) println("NaN Acceptance $acceptanceNumerator / $acceptanceDenominator logPiv = ${state.logProbOfPivotState} transition prob = ${transitionProb(revertState.reversePivot)} columnWeight = ${columnWeights.P(revertState.reversePivot.col)} nPivots = ${nPivots(revertState.reversePivot.col)}");
    std::cout << "Log acceptance is " << destinationProb << " - " << sourceProb << " + " << proposalPivot.logAcceptanceContribution << " = " << logAcceptance << std::endl;
    if (std::isnan(logAcceptance) || Random::nextDouble() <= exp(std::min(0.0,logAcceptance))) { // explicity accept if both numerator and denominator are -infinity
        // Accept proposal
        std::cout << "Accepting" << std::endl;
        pivot(proposalPivot);
    } else {
        // reject
        std::cout << "Rejecting" << std::endl;
    }
}

//double SimplexMCMC::logTransitionProb(const ColumnPivot &proposal) {
//    return -log(nNonBasic())-log(proposal.pivotRows.size());
//}
//
//double SimplexMCMC::logReverseTransitionProb(const ColumnPivot &proposal) {
//    int degeneracyCount = 1;
//    for(int i=1; i<proposal.col.size(); ++i) {
//        if(i != proposal.i && fabs(proposal.col[i]) > 1e-6 && (b[i] == u[head[i]] || b[i] == l[head[i]])) ++degeneracyCount;
//    }
//    return -log(nNonBasic()) - log(degeneracyCount);
//}

void SimplexMCMC::randomWalk() {
    pivot(proposePivot());
}

ProbabilisticColumnPivot SimplexMCMC::proposePivot() {
//    return ProbabilisticColumnPivot(*this, proposeColumn());
    return ProbabilisticColumnPivot(*this);
}

// To be based on rate of change of L1-norm infeasibility objective?
// uniform prob for now
int SimplexMCMC::proposeColumn() {
    // choose a pivot column
    return Random::nextInt(1,n - m + 1);
}


// Moved to ProbabilisticColumnPivot
//void SimplexMCMC::proposeRow(ColumnPivot &colProposal) {
//    if(colProposal.pivotRows.size() > 0) {
//        colProposal.i = colProposal.pivotRows[Random::nextInt(colProposal.pivotRows.size())];
//    } else {
//        colProposal.i = -1; // incoming var goes to opposite limit.
//    }
//}


// pivots-in any auxiliary variables if there exists a basic variable on a bound
// that has non-zero value in the incoming column.
// Uses a greedy algorithm, pivoting in the candidate variable
void SimplexMCMC::toCanonicalState() {
    std::vector<int> nonBasicAux = auxiliaries();
    bool pivoted = false;
    do {

    } while(pivoted);
}

// returns the set of nonBasic variables (by j-value) that correspond
// to auxiliary variables in the original problem.
std::vector<int> SimplexMCMC::auxiliaries() {
    std::vector<int> nbAux;
    for(int j = 1; j<=nNonBasic(); ++j) {
        if(kSimTokProb[head[j]] <= originalProblem.nConstraints()) nbAux.push_back(j);
    }
    return nbAux;
}

int SimplexMCMC::countFractionalPivCols() {
    int nFractionals = 0;
    double delta;
    for(int j=1; j<=nNonBasic(); ++j) {
        delta = fabs(ColumnPivot(*this, j).deltaj);
        if(delta < 0.9999 && delta > 0.0001) ++nFractionals;
    }
    return nFractionals;
}

void SimplexMCMC::updateLPSolution(const ProposalPivot &pivot) {
    int nConstraints = originalProblem.nConstraints();
    for(int i=1; i<=nBasic(); ++i) {
        int kProb = kSimTokProb[head[i]];
        if(kProb > nConstraints) lpSolution[kProb - nConstraints] += pivot.col[i] * pivot.deltaj;
    }
    int kCol = kSimTokProb[head[nBasic() + pivot.j]];
    if(kCol > nConstraints) lpSolution[kCol - nConstraints] += pivot.deltaj;
}

void SimplexMCMC::revertLPSolution(const ProposalPivot &pivot) {
    int nConstraints = originalProblem.nConstraints();
    int kCol = kSimTokProb[head[nBasic() + pivot.j]];
    if(kCol > nConstraints) lpSolution[kCol - nConstraints] -= pivot.deltaj;
    for(int i=1; i<=nBasic(); ++i) {
        int kProb = kSimTokProb[head[i]];
        if(kProb > nConstraints) lpSolution[kProb - nConstraints] -= pivot.col[i] * pivot.deltaj;
    }
}

bool SimplexMCMC::solutionIsPrimaryFeasible() {
    const double tol = 1e-8;
    for(int i=1; i<nBasic(); ++i) {
        int k = head[i];
        if(b[i] < l[k] - tol || b[i] > u[k] + tol) return false;
    }
    return true;
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