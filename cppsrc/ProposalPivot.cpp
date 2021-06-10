//
// Created by daniel on 13/05/2021.
//

#include <cmath>
#include "ProposalPivot.h"

//ProposalPivot::ProposalPivot(glp::Simplex &lp, int j): j(j), col(lp.tableauCol(j)) {
//    int kIncoming = lp.head[lp.m + j];
//    double deltaMin = lp.u[kIncoming] - lp.l[kIncoming];
//    double DXi;
//    double DXj;
//    for(int i=1; i<col.size();++i) {
//        if(col[i] != 0.0) {
//            int kOutgoing = lp.head[i];
//            bool outgoingToUpperBound = (col[i] > 0.0) ^ lp.isAtUpperBound(j);
//            if(outgoingToUpperBound) {
//                DXi = lp.u[kOutgoing] - lp.b[i];
//            } else {
//                DXi = lp.l[kOutgoing] - lp.b[i];
//            }
//            DXj = fabs(DXi/col[i]);
//            if (DXj < deltaMin) {
//                pivotRows.clear();
//                pivotRows.push_back(i);
//                deltaMin = DXj;
//            } else if (DXj == deltaMin) {
//                pivotRows.push_back(i);
//            }
//        }
//    }
//    if(pivotRows.empty()) this->i = -1; else this->i = pivotRows.front();
//}


// returns the column for the reverse pivot
// or the empty vector if i < 1 (i.e. if the pivot is just a bound swap of the j'th column)
std::vector<double> ProposalPivot::reverseCol() const {
    std::vector<double> revCol(col.size());
    if(i < 1) {
        for(int irev=0; irev < col.size(); ++irev) {
            revCol[irev] = col[irev];
        }
    } else {
        double denominator = col[i];
        for(int irev=0; irev < col.size(); ++irev) {
            revCol[irev] = col[irev] / denominator;
        }
        revCol[i] /= denominator;
    }
    return revCol;
}
