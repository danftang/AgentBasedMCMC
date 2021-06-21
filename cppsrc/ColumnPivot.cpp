//
// Created by daniel on 17/05/2021.
//

#include <cmath>
#include <cassert>
#include "ColumnPivot.h"
#include "StlStream.h"

ColumnPivot::ColumnPivot(glp::Simplex &lp, int j, std::vector<double> column): ProposalPivot(-1, j, std::move(column)) {
    int kIncoming = lp.head[lp.m + j];
    deltaj = lp.u[kIncoming] - lp.l[kIncoming];
    double DXi;
    double DXj;
    for(int i=1; i<col.size();++i) { // TODO: allow optionally degenerate pivots if both incoming and outgoing are on bounds?
        if(fabs(col[i]) > tol) {
            int kOutgoing = lp.head[i];
            bool outgoingToUpperBound = (col[i] > 0.0) ^ lp.isAtUpperBound(j);
            if(outgoingToUpperBound) {
                DXi = lp.u[kOutgoing] - lp.b[i];
            } else {
                DXi = lp.l[kOutgoing] - lp.b[i];
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
    if(lp.isAtUpperBound(j)) deltaj = -deltaj;
    orderPivotRows(lp);
}


ColumnPivot ColumnPivot::reverse(glp::Simplex &lp) const {
    ColumnPivot revCol(lp, j, reverseCol());
    revCol.i = i;
    return revCol;
}


// orders pivot rows so that all structural vars in the original LP come before auxiliary vars
void ColumnPivot::orderPivotRows(glp::Simplex &lp) {
    nStructuralPivotRows = 0;
    while(lp.kSimTokProb[lp.head[pivotRows[nStructuralPivotRows]]] > lp.originalProblem.nConstraints() && nStructuralPivotRows < pivotRows.size()) {
        ++nStructuralPivotRows;
    }
    for(int entry = nStructuralPivotRows + 1; entry < pivotRows.size(); ++entry) {
        int i = pivotRows[entry];
        if(lp.kSimTokProb[lp.head[i]] > lp.originalProblem.nConstraints()) {
            // swap current entry with first non-structural
            pivotRows[entry] = pivotRows[nStructuralPivotRows];
            pivotRows[nStructuralPivotRows++] = i;
        }
    }
}

