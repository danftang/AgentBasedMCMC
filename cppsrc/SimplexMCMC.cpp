//
// Created by daniel on 15/04/2021.
//

#include "SimplexMCMC.h"
#include "ProposalPivot.h"
#include "Random.h"
#include "Phase2Pivot.h"
#include "PotentialEnergyPivot.h"
#include "ConvexPMF.h"
#include "debug.h"
#include <algorithm>
#include <cmath>
#include <cfloat>


// Initial basis of 'prob' should contain no fixed vars and all remaining auxiliary vars
// use
SimplexMCMC::SimplexMCMC(
        const glp::Problem &prob,
        std::function<double (const std::vector<double> &)> logProb,
        const std::vector<double> &initialState)
        : Simplex(prob),
        logProbFunc(std::move(logProb)),
        proposalFunction(*this) {
    setObjective(glp::SparseVec());
    if(initialState.size() != 0) setLPState(initialState);
    findFeasibleStartPoint();
    proposalFunction.initCache();
}


// Starting with zero exactEndState, find a exactEndState that satisfies the observations and
// the interaction rules.
void SimplexMCMC::findFeasibleStartPoint() {
    int iterations = 0;
    while(!solutionIsPrimaryFeasible()) {
        Phase1Pivot proposal(*this);
        ProposalPivot & p = proposal;
        pivot(proposal);
//        std::cout << iterations << " Pivoted on " << proposal.i << ", " << proposal.j << " " << proposal.deltaj << " " << proposal.leavingVarToUpperBound
//                  << " " << isAtUpperBound(proposal.j)
//                  << "  Infeasibility = " << infeasibility() << std::endl;
        debug(if(iterations%256 == 0) std::cout << "phase1 iteration:" << iterations << " infeasibility = " << infeasibility() << std::endl);
        iterations++;

    } ;
    debug(std::cout << "Found initial exactEndState in " << iterations << " iterations" << std::endl);
}


// Calculate the degeneracy probability of the current state
// from a cold start (i.e. nothing already calculated)
//
// If we split the variables into the ones that are on a bound, d, and those that aren't, x,
// then any basis that includes all x variables will have the same exactEndState.
// However, on an integer exactEndState, all variables are on their bound, so any
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

//double SimplexMCMC::lnFractionalPenalty() {
//    // fractions must be in the basis, so check b
//    double penalty = 0.0;
//    for(int i=1; i<=m; ++i) {
//        if(fabs(round(beta[i]) - beta[i]) > tol) {
//            penalty += fractionalK;
//        }
//    }
//    return penalty;
//}


// Choose a pivot
// and reject based on Metropolis-Hastings
// returns the next nextSample
const std::vector<double> & SimplexMCMC::nextSample() {
    int infeasibleCount = 0;
    bool sampleIsFeasible;
    do {
        const ProposalPivot &proposalPivot = proposePivot();
        lastSampleWasAccepted = processProposal(proposalPivot);
        sampleIsFeasible = solutionIsPrimaryFeasible();
        if(sampleIsFeasible) {
            feasibleStatistics.update(lastSampleWasAccepted, proposalPivot);
        } else {
            infeasibleStatistics.update(lastSampleWasAccepted, proposalPivot);
            debug(if(infeasibleCount%1000 == 100) std::cout << infeasibleCount << " Infeasibility = " << infeasibility() << std::endl);
            ++infeasibleCount;
        }
    } while(!sampleIsFeasible);
    return X();
}


bool SimplexMCMC::processProposal(const ProposalPivot &proposalPivot) {
//    std::cout << "Processing proposal " << proposalPivot.i << ", " << proposalPivot.j
//              << " leavingVarToUpperBound = " << proposalPivot.leavingVarToUpperBound
//              << " deltaj = " << proposalPivot.deltaj
//              << " incoming var = " << (isAtUpperBound(proposalPivot.j)?u[head[nBasic()+proposalPivot.j]]:l[head[nBasic()+proposalPivot.j]]);
//    if(proposalPivot.i > 0) {
//        std::cout << " b[i] = " << b[proposalPivot.i] << " leaving var limits " << l[head[proposalPivot.i]]
//                  << ":" << u[head[proposalPivot.i]];
//    }
//    std::cout << std::endl;

    if(proposalPivot.i > 0  && kSimTokProb[head[proposalPivot.i]] <= nBasic()) std::cout << "Warning: Proposing auxiliary var to leave basis" << std::endl;

    double sourceProb = logProbFunc(X());
//    std::cout << "Source LP state is: " << glp::SparseVec(lpSolution) << std::endl;
    updateLPSolution(proposalPivot);
    double destinationProb = logProbFunc(lpSolution);
//    std::cout << "Destination LP state is: " << glp::SparseVec(lpSolution) << std::endl;
    double logAcceptance = destinationProb - sourceProb + proposalPivot.logAcceptanceContribution;
//    debug(
//            if(std::isnan(logAcceptance)) std::cout << "Log acceptance is " << logAcceptance << std::endl
//            );
//    if(isnan( logAcceptance )) println("NaN Acceptance $acceptanceNumerator / $acceptanceDenominator logPiv = ${state.logProbOfPivotState} transition prob = ${transitionProb(revertState.reversePivot)} columnWeight = ${columnWeights.P(revertState.reversePivot.col)} nPivots = ${nPivots(revertState.reversePivot.col)}");
//    debug(std::cout << "Log acceptance is " << destinationProb << " - " << sourceProb << " + " << proposalPivot.logAcceptanceContribution << " = " << logAcceptance << std::endl);
    if (std::isnan(logAcceptance) || Random::nextDouble() <= exp(std::min(0.0,logAcceptance))) { // explicity accept if both numerator and denominator are -infinity
        // Accept proposal
//        debug(std::cout << "Accepted proposal with log acceptance = " << logAcceptance << std::endl);
        //        std::cout << "Accepting deltaj = " << proposalPivot.deltaj << std::endl;
        assert(proposalPivot.i == -1 || isStructural(head[proposalPivot.i]));
        assert(proposalPivot.i == -1 || proposalPivot.col[proposalPivot.i] > 0.9 || proposalPivot.col[proposalPivot.i] < -0.9);
//        if(proposalPivot.deltaj != 0.0) {
//            if(proposalPivot.i > 0) {
//                std::cout << "Pivoting " << proposalPivot.deltaj << " " << proposalPivot.col[proposalPivot.i] << " "
//                          << proposalPivot.i << " / " << nBasic() << " " << proposalPivot.j << " / " << nNonBasic()
//                          << std::endl;
//            } else {
//                std::cout << "Swapping " << proposalPivot.deltaj << " "
//                          << proposalPivot.j << " / " << nNonBasic()
//                          << std::endl;
//
//            }
//        }
        pivot(proposalPivot);
        return true; // accept
    }
//    debug(std::cout << "Rejecting" << std::endl);
    revertLPSolution(proposalPivot);
    return false; // reject
}

//double SimplexMCMC::logTransitionProb(const Phase2Pivot &proposal) {
//    return -log(nNonBasic())-log(proposal.pivotRows.size());
//}
//
//double SimplexMCMC::logReverseTransitionProb(const Phase2Pivot &proposal) {
//    int degeneracyCount = 1;
//    for(int i=1; i<proposal.col.size(); ++i) {
//        if(i != proposal.i && fabs(proposal.col[i]) > 1e-6 && (b[i] == u[head[i]] || b[i] == l[head[i]])) ++degeneracyCount;
//    }
//    return -log(nNonBasic()) - log(degeneracyCount);
//}


// lpState should contain the values of the strucutral variables of the original LP
void SimplexMCMC::setLPState(const std::vector<double> &lpState) {
    for(int j=1; j <= nNonBasic(); ++j) {
        int kSim = head[nBasic() + j];
        assert(isStructural(kSim));
        int kLP = kSimTokProb[kSim];
        double v = lpState[kLP - nBasic()];
        flag[j] = (fabs(u[kSim] - v) < fabs(l[kSim] - v));
    }
    spx_eval_beta(this, beta);
}

void SimplexMCMC::randomWalk() {
    pivot(Phase2Pivot(*this, Random::nextInt(1, nNonBasic()+1)));
}

const ProposalPivot &SimplexMCMC::proposePivot() {
    return proposalFunction.nextProposal();
}

// To be based on rate of change of L1-norm infeasibility objective?
// uniform prob for now
//int SimplexMCMC::proposeColumn() {
//    // nextIntFromDiscrete a pivot column
//    return Random::nextInt(1,n - m + 1);
//}


// Moved to ProbabilisticColumnPivot
//void SimplexMCMC::proposeRow(Phase2Pivot &colProposal) {
//    if(colProposal.pivotRows.size() > 0) {
//        colProposal.i = colProposal.pivotRows[Random::nextInt(colProposal.pivotRows.size())];
//    } else {
//        colProposal.i = -1; // incoming var goes to opposite limit.
//    }
//}


// pivots-in any auxiliary variables if there exists a basic variable on a bound
// that has non-zero value in the incoming column.
// Uses a greedy algorithm, pivoting in the candidate variable
//void SimplexMCMC::toCanonicalState() {
//    std::vector<int> nonBasicAux = auxiliaries();
//    bool pivoted = false;
//    do {
//
//    } while(pivoted);
//}


// returns the set of nonBasic variables (by j-value) that correspond
// to auxiliary variables in the original problem.
//std::vector<int> SimplexMCMC::auxiliaries() {
//    std::vector<int> nbAux;
//    for(int j = 1; j<=nNonBasic(); ++j) {
//        if(kSimTokProb[head[nBasic() + j]] <= originalProblem.nConstraints()) nbAux.push_back(j);
//    }
//    return nbAux;
//}


int SimplexMCMC::countFractionalPivCols() {
    int nFractionals = 0;
    double delta;
    for(int j=1; j<=nNonBasic(); ++j) {
        delta = fabs(Phase2Pivot(*this, j).deltaj);
        if(delta < 0.9999 && delta > 0.0001) ++nFractionals;
    }
    return nFractionals;
}


int SimplexMCMC::infeasibilityCount() {
    int infeasibility = 0;
    for(int i=1; i<=nBasic();++i) {
        int k = head[i];
        if(beta[i] < l[k] - tol || beta[i] > u[k] + tol) ++infeasibility;
    }
    return infeasibility;
}


double SimplexMCMC::infeasibility() {
    double dist = 0.0;
    for(int i=1; i <= nBasic(); ++i) {
        int k = head[i];
        double v = beta[i];
        if(v < l[k]) {
//            std::cout << "l[k]= " << l[k] << " v = " << v << std::endl;
            dist += l[k] - v;
        } else if(v > u[k]) {
//            std::cout << "v =  " << v << " u[k] =  " << u[k] << std::endl;
            dist += v - u[k];
        }
    }
    return dist;
}


void SimplexMCMC::updateLPSolution(const ProposalPivot &pivot) {
    int nConstraints = nBasic();
    for(int i=1; i<=nBasic(); ++i) {
        int kProb = kSimTokProb[head[i]];
        if(kProb > nConstraints) lpSolution[kProb - nConstraints] += pivot.col[i] * pivot.deltaj;
    }
    int kCol = kSimTokProb[head[nBasic() + pivot.j]];
    if(kCol > nConstraints) lpSolution[kCol - nConstraints] += pivot.deltaj;
    lpSolutionIsValid(false);
}


void SimplexMCMC::revertLPSolution(const ProposalPivot &pivot) {
    int nConstraints = nBasic();
    int kCol = kSimTokProb[head[nBasic() + pivot.j]];
    if(kCol > nConstraints) lpSolution[kCol - nConstraints] -= pivot.deltaj;
    for(int i=1; i<=nBasic(); ++i) {
        int kProb = kSimTokProb[head[i]];
        if(kProb > nConstraints) lpSolution[kProb - nConstraints] -= pivot.col[i] * pivot.deltaj;
    }
    lpSolutionIsValid(true);
}


bool SimplexMCMC::solutionIsPrimaryFeasible() {
    for(int i=1; i<=nBasic(); ++i) {
        int k = head[i];
        if(beta[i] < (l[k] - tol) || beta[i] > (u[k] + tol)) return false;
    }
    return true;
}

bool SimplexMCMC::solutionIsInteger() {
    for(int i=1; i<=nBasic(); ++i) {
        if(fabs(beta[i] - round(beta[i])) > tol) return false;
    }
    return true;
};

glp::Problem &SimplexMCMC::initialiseProblem(glp::Problem &lp) {
    lp.advBasis();
    lp.warmUp();
    return(lp);
}


void SimplexMCMC::SampleStatistics::update(bool accepted, const ProposalPivot &proposal) {
    ++nSamples;
    if(accepted) {
        ++nAccepted;
        if(fabs(proposal.deltaj) < tol) {
            if(proposal.i < 1) ++nNulls;
        } else {
            ++nNonDegenerate;
            if(proposal.i < 1) ++nSwaps;
        }
    }
}


// Check initial basis contains no fixed vars and all auxiliaries are in the basis
bool SimplexMCMC::abmSanityChecks() {
    for(int k=1; k<=nVars(); ++k) {
        if(l[k] == u[k] && isStructural(k)) {
            std::cerr << "WARNING: Non-auxiliary, fixed variable in the Simplex" << k << std::endl;
            return false;
        }
        if(isStructural(k) && (l[k] != 0.0 || u[k] != 1.0 )) {
            std::cerr << "WARNING: non-binary structural var " << k << std::endl;
            return false;
        }
    }
    for(int j=1; j<=nNonBasic(); ++j) {
        int k = head[nBasic()+j];
        if(!isStructural(k)) {
            std::cerr << "WARNING: Non-basic auxiliary variable " << j << std::endl;
            return false;
        }
        if(l[k] == -DBL_MAX || u[k] == DBL_MAX) {
            std::cerr << "WARNING: Non-basic variable with infinite bound " << j << std::endl;
            return false;
        }
    }
//    assert(originalProblem.isValidSolution(X()));
    return true;
}


std::ostream &operator <<(std::ostream &out, const SimplexMCMC::SampleStatistics &stats) {
    out << "Total samples           " << stats.nSamples << std::endl;
    out << "accepted/total          " << stats.nAccepted*100.0/stats.nSamples << "%" << std::endl;
    out << "non-degenerate/accepted " << stats.nNonDegenerate*100.0/stats.nAccepted << "%" << std::endl;
    out << "pivots/non-degenerate   " << (stats.nNonDegenerate-stats.nSwaps)*100.0/stats.nNonDegenerate << "%" << std::endl;
    out << "swaps/non-degenerate    " << stats.nSwaps*100.0/stats.nNonDegenerate << "%" << std::endl;
    out << "nulls/degenerate        " << stats.nNulls*100.0/(stats.nAccepted-stats.nNonDegenerate) << "%" << std::endl;
    return out;
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