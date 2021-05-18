//
// Created by daniel on 17/05/2021.
//

#include <cmath>
#include <cassert>
#include "ColumnPivot.h"
#include "StlStream.h"

ColumnPivot::ColumnPivot(glp::Simplex &lp, int j, std::vector<double> column): Pivot(-1,j,column) {
    int kIncoming = lp.head[lp.m + j];
    double deltaMin = lp.u[kIncoming] - lp.l[kIncoming];
    double DXi;
    double DXj;
    for(int i=1; i<col.size();++i) {
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
    if(!pivotRows.empty()) this->i = pivotRows.front();
}


ColumnPivot ColumnPivot::reverse(glp::Simplex &lp) const {
    ColumnPivot revCol(lp, j, reverseCol());
    revCol.i = i;
    return revCol;
}
