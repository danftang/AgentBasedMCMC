//
// Created by daniel on 28/06/2021.
//

#include "Phase1Pivot.h"
#include <float.h>
#include <cassert>

Phase1Pivot::Phase1Pivot(glp::Simplex &simplex): ProposalPivot(simplex) {
    chooseCol();
    chooseRow();
}


void Phase1Pivot::chooseCol() {
    setSimplexToInfeasibilityObjective();
    simplex.recalculatePi();
    setCol(*std::min_element(simplex.pi.begin(), simplex.pi.end()));
}


void Phase1Pivot::chooseRow() {
    std::multimap<double,int> pivots = getPivotsByDeltaJ();

    auto pivotIt = pivots.begin();
    int pivotIndex = pivotIt->second;
    double dDf_dDj = colInfeasibilityGradient(pivots.begin()->first - tol);
    while( dDf_dDj < 0.0) {
        ++pivotIt;
        pivotIndex = pivotIt->second;
        if(pivotIndex < nonZeroRows.size() * 2) {
            dDf_dDj += fabs(col[nonZeroRows[pivotIndex / 2]]);
        } else {
            dDf_dDj += 1.0;
        }
    }
    setToPivotIndex(pivotIndex);
}


void Phase1Pivot::setToPivotIndex(int pivotIndex) {
    leavingVarToUpperBound = pivotIndex % 2;
    if(pivotIndex < 2 * nonZeroRows.size()) {
        i = nonZeroRows[pivotIndex / 2];
        int leavingk = simplex.head[i];
        deltaj = ((leavingVarToUpperBound ? simplex.u[leavingk] : simplex.l[leavingk]) - simplex.b[i]) / col[i];
    } else {
        i = -1; // column does bound swap.
        int k = simplex.head[simplex.nBasic() + j];
        deltaj = leavingVarToUpperBound?(simplex.u[k]-simplex.l[k]):(simplex.l[k]-simplex.u[k]);
    }
}


int Phase1Pivot::setSimplexToInfeasibilityObjective() {
    // set objective to out-of-bounds rows
    int infeasibilityCount = 0;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        if(simplex.b[i] < simplex.l[k] - tol) {
            simplex.c[i] = -1.0;
            ++infeasibilityCount;
        } else if(simplex.b[i] > simplex.u[k] + tol) {
            simplex.c[i] = 1.0;
            ++infeasibilityCount;
        } else {
            simplex.c[i] = 0.0;
        }
    }
    return infeasibilityCount;
}


// returns gradient after perturbation of this column by deltaj
double Phase1Pivot::colInfeasibilityGradient(double deltaj) {
    double grad = 0.0;
    for(int nzi : nonZeroRows) {
        int rowk = simplex.head[nzi];
        grad += col[nzi] * infeasibilityGradient(simplex.b[nzi] + col[nzi] * deltaj,
                                                 simplex.l[rowk],
                                                 simplex.u[rowk]);
    }
    int colk = simplex.head[simplex.nBasic() + j];
    double xUpperBound = simplex.u[colk];
    double xLowerBound = simplex.l[colk];
    double xj = simplex.isAtUpperBound(j)?xUpperBound:xLowerBound;
    grad += infeasibilityGradient(xj + deltaj, xLowerBound, xUpperBound);
    return grad;
}


double Phase1Pivot::infeasibilityGradient(double v, double lowerBound, double upperBound) {
    if(v >= upperBound) {
        return 1.0;
    } else if(v < lowerBound) {
        return -1.0;
    }
    return 0.0;
}
