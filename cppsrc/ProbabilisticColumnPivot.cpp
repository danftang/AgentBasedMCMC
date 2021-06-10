//
// Created by daniel on 08/06/2021.
//

#include <cmath>
#include "ProbabilisticColumnPivot.h"
#include "Random.h"

// set up cumulative probability
ProbabilisticColumnPivot::ProbabilisticColumnPivot(glp::Simplex &simplex, int j, std::vector<double> column): ProposalPivot(-1, j, std::move(column)) {

    std::multimap<double, int> transitions; // from delta_j to PMF-index.

    // identify active rows (non-zero in structural row)
    for(int i=1; i<col.size(); ++i) {
        int k = simplex.head[i];
        if(fabs(col[i]) > tol && simplex.kSimTokProb[k] > simplex.originalProblem.nConstraints()) activeRows.push_back(i);
    }

    // calculate delta_js for each possible pivot
    for(int m=0; m<activeRows.size(); ++m) {
        int i = activeRows[m];
        int k = simplex.head[i];
        double deltaLB = (simplex.b[i] - simplex.l[k])/col[i];
        transitions.emplace(deltaLB, 2*m);
        double deltaUB = (simplex.b[i] - simplex.u[k])/col[i];
        transitions.emplace(deltaUB, 2*m + 1);
    }

    // add entry for bound swap of column
    int colk = simplex.head[j];
    double boundSwapDelta = simplex.isAtUpperBound(j)?(simplex.l[colk]- simplex.u[colk]):(simplex.u[colk]- simplex.l[colk]);
    transitions.emplace(boundSwapDelta, 2*activeRows.size());

    // now populate pivotPMF
    pivotPMF.resize(activeRows.size() * 2 + 1, 0.0);
    auto firstPositivePivot = transitions.lower_bound(0.0);
    // first do +ve Delta_j pivots (if any)
    double DeltaF = 0.0;
    double lastDj = 0.0;
    double dDf_dDj = colFeasibilityGradient(simplex, true);
    for(auto piv = firstPositivePivot; piv != transitions.end(); ++piv) {
        auto [Dj, pmfIndex] = *piv;
        DeltaF += (Dj-lastDj)*dDf_dDj;
        lastDj = Dj;
        dDf_dDj += 1.0;
        // only pivot on unity pivot points (this ensures solutions remain integer if coeffs are integer to start with)
        // TODO: Understand consequences of this design decision
        pivotPMF[pmfIndex] = (pmfIndex == 2*activeRows.size() || fabs(fabs(col[activeRows[pmfIndex/2]])-1.0) < tol )?exp(kappa * DeltaF):0.0;
    }
    // now do -ve Delta_j pivots
    DeltaF = 0.0;
    lastDj = 0.0;
    dDf_dDj = colFeasibilityGradient(simplex, false);
    for(auto piv = std::make_reverse_iterator(firstPositivePivot); piv != transitions.rend(); ++piv) {
        auto [Dj, pmfIndex] = *piv;
        DeltaF += (Dj-lastDj)*dDf_dDj;
        lastDj = Dj;
        dDf_dDj -= 1.0;
        // only pivot on unity pivot points (this ensures solutions remain integer if coeffs are integer to start with)
        pivotPMF[pmfIndex] = (pmfIndex == 2*activeRows.size() || fabs(fabs(col[activeRows[pmfIndex/2]])-1.0) < tol )?exp(kappa * DeltaF):0.0;
    }

    // choose row
    int pivotChoice = Random::choose(pivotPMF.begin(), pivotPMF.end());
    if(pivotChoice < 2*activeRows.size()) {
        i = activeRows[pivotChoice / 2];
        leavingVarToUpperBound = pivotChoice % 2;
        int leavingk = simplex.head[i];
        delta = ((leavingVarToUpperBound?simplex.u[leavingk]:simplex.l[leavingk])  - simplex.b[i])/col[i];
    } else {
        i = -1; // column does bound swap.
        leavingVarToUpperBound = !simplex.isAtUpperBound(j);
        delta = boundSwapDelta;
    }
    logTransitionRatio = 0.0;
}


// returns the gradient of the feasibility function with respect to changes in the value of
// this column in either the forward (increasing val) or backward (decreasing val) directions.
// (since the col is on a boundary, the gradient is discontinuous at this point so direction
// must be specified).
double ProbabilisticColumnPivot::colFeasibilityGradient(glp::Simplex &lp, bool forward) {
    double grad = 0.0;
    for(int nzi : activeRows) {
        grad += iFeasibilityGradient(lp, nzi, forward);
    }
    if(lp.isAtUpperBound(j)) {
        if(forward) grad += 1.0;
    } else {
        if(!forward) grad -= 1.0;
    }
    return grad;
}

// returns the gradient of the feasibility function for the current value of the i'th row.
// If 'forward' is true, returns the gradient in the forward direction, otherwise returns the grad
// in the backward direction (these can be different at boundaries since the gradient
// is discontinuous there).
double ProbabilisticColumnPivot::iFeasibilityGradient(glp::Simplex &lp, int i, bool forward) {
    int k = lp.head[i];
    double v = lp.b[i];
    if(double upperBound = lp.u[k]; v > upperBound - forward?tol:-tol) {
        return 1.0;
    } else if(double lowerBound = lp.l[k]; v < lowerBound - forward?tol:-tol) {
        return -1.0;
    }
    return 0.0;
}
