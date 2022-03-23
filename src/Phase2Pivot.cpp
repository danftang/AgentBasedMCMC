//
// Created by daniel on 17/05/2021.
//

#include <cmath>
#include <cassert>
#include "Phase2Pivot.h"
#include "StlStream.h"
#include "SimplexMCMC.h"

Phase2Pivot::Phase2Pivot(SimplexMCMC &simplex, int j) : Phase2Pivot(simplex, j, simplex.tableauCol(j)) { }

Phase2Pivot::Phase2Pivot(SimplexMCMC &lp, int j, std::vector<double> column): ProposalPivot(lp, -1, j, std::move(column)) {
    int kIncoming = lp.head[lp.m + j];
    deltaj = lp.u[kIncoming] - lp.l[kIncoming];
    double DXi;
    double DXj;
    for(int i=1; i<col.vec.size();++i) {
        if(fabs(col[i]) > tol) {
            int kOutgoing = lp.head[i];
            bool outgoingToUpperBound = (col[i] > 0.0) ^ lp.isAtUpperBound(j);
            if(outgoingToUpperBound) {
                DXi = lp.u[kOutgoing] - lp.beta[i];
            } else {
                DXi = lp.l[kOutgoing] - lp.beta[i];
            }
            DXj = fabs(DXi/col[i]);
            if (DXj < deltaj) {
                pivotRows.clear();
                pivotRows.push_back(i);
                deltaj = DXj;
            } else if (DXj == deltaj) {
                pivotRows.push_back(i);
            }
        }
    }
    if(!pivotRows.empty()) this->i = pivotRows.front();
    if(lp.isAtUpperBound(j)) {
        deltaj = -deltaj;
        leavingVarToUpperBound = (col[i] < 0.0);
    } else {
        leavingVarToUpperBound = (col[i] > 0.0);
    }
//    orderPivotRows(lp);
}


//Phase2Pivot Phase2Pivot::reverse(SimplexMCMC &lp) const {
//    Phase2Pivot revCol(lp, j, reverseCol());
//    revCol.i = i;
//    return revCol;
//}


// orders pivot rows so that all structural vars in the original LP come before auxiliary vars
void Phase2Pivot::orderPivotRows(SimplexMCMC &lp) {
    nStructuralPivotRows = 0;
    while(nStructuralPivotRows < pivotRows.size() && lp.isStructural(lp.head[pivotRows[nStructuralPivotRows]])) {
        ++nStructuralPivotRows;
    }
    for(int entry = nStructuralPivotRows + 1; entry < pivotRows.size(); ++entry) {
        int i = pivotRows[entry];
        if(lp.isStructural(lp.head[i])) {
            // swap current entry with first non-structural
            pivotRows[entry] = pivotRows[nStructuralPivotRows];
            pivotRows[nStructuralPivotRows++] = i;
        }
    }
}

