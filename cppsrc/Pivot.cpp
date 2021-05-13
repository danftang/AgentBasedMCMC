//
// Created by daniel on 13/05/2021.
//

#include <cmath>
#include "Pivot.h"

Pivot::Pivot(glp::Simplex &lp, int j): Pivot(lp, 0, j, lp.tableauCol(j)) { }

Pivot::Pivot(glp::Simplex &lp, int i, int j, std::vector<double> column): i(i), j(j), col(std::move(column)) {
    int kIncoming = lp.head[lp.m + j];
    double deltaMin = lp.u[kIncoming] - lp.l[kIncoming];
    double DXi;
    double DXj;
    for(int i=1; i<=lp.nRows();++i) {
        if(col[i] != 0.0) {
            int kOutgoing = lp.head[i];
            bool outgoingToUpperBound = (col[i] > 0.0) ^ lp.isAtUpperBound(j);
            if(outgoingToUpperBound) {
                DXi = lp.u[kOutgoing] - lp.b[i];
            } else {
                DXi = lp.l[kOutgoing] - lp.b[i];
            }
            DXj = fabs(DXi/col[i]);
            if (DXj < deltaMin) {
                pivotRows.clear();
                pivotRows.push_back(i);
                deltaMin = DXj;
            } else if (DXj == deltaMin) {
                pivotRows.push_back(i);
            }
        }
    }
    if(i<1 && !pivotRows.empty()) {
        this->i = pivotRows.front();
    }
}

std::vector<double> Pivot::reverseCol() {
    std::vector<double> revCol;
    revCol.resize(col.size());
    double denominator = fabs(col[i]);
    for(int i=0; i<col.size(); ++i) {
        revCol[i] = col[i] / denominator;
    }
    return revCol;
}
