//
// Created by daniel on 13/05/2021.
//

#include <cmath>
#include <cassert>
#include <float.h>
#include "ProposalPivot.h"
#include "SimplexMCMC.h"


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

void ProposalPivot::initNonZeroRows() {
    nonZeroRows.clear();
    for(int i=1; i<col.size(); ++i) {
        if(fabs(col[i]) > tol) nonZeroRows.push_back(i);
    }
}


// calculate a sorted map from delta_js to pivot-index for each possible pivot,
// where pivot index is 2*nonZeroRowsIndex + leavingVarToUpperBound
// if nonZeroRowsIndex == nonZeroRows.size() then the pivot is a bound swap
// on this column, the basis is unchanged and leavingVarToUpperBound
// refers to this column. Note that the null pivot (this column stays on
// its current bound) is included in order to help calculation of gradients.
std::multimap<double,int> ProposalPivot::getPivotsByDeltaJ() {
    assert(j > 0);

    std::multimap<double, int> allPivots; // from delta_j to LogPMF-index.

    for(int m=0; m < nonZeroRows.size(); ++m) {
        int i = nonZeroRows[m];
        int k = simplex.head[i];
        double rowLowerBound = simplex.l[k];
        double rowUpperBound = simplex.u[k];
        assert(col[i] != 0.0);
        if(rowLowerBound > -DBL_MAX) {
            double deltaLB = (rowLowerBound - simplex.beta[i]) / col[i];
            allPivots.emplace(deltaLB, 2 * m);
        }
        if(rowUpperBound < DBL_MAX) {
            double deltaUB = (rowUpperBound - simplex.beta[i]) / col[i];
            allPivots.emplace(deltaUB, 2 * m + 1);
        }
    }

    // add entries for upper and lower bounds of this column
    int k = simplex.head[simplex.nBasic() + j];
    double Xj = simplex.isAtUpperBound(j)?simplex.u[k]:simplex.l[k];
    if(simplex.l[k] > -DBL_MAX) allPivots.emplace(simplex.l[k] - Xj, 2 * nonZeroRows.size());
    if(simplex.u[k] < DBL_MAX) allPivots.emplace(simplex.u[k] - Xj, 2 * nonZeroRows.size() + 1);

    return allPivots;
}


std::multimap<double,int> ProposalPivot::getPivotsByInfeasibility() {
    std::multimap<double,int> pivotsByDeltaj = getPivotsByDeltaJ();
    std::multimap<double,int> pivotsByInf;
//    assert(simplex.reducedObjective(j) == colInfeasibilityGradient(0.0)); // TEST

    auto pivotIt = pivotsByDeltaj.begin();
    double lastDj = pivotIt->first;
    int minPivotIndex = pivotIt->second;
    double infeas = 0.0; // infeasibility(lastDj);
    double minDeltaF = infeas;
    double dDf_dDj = colInfeasibilityGradient(lastDj - tol);
    for(auto [Dj, pivotId] : pivotsByDeltaj) {
        infeas += (Dj - lastDj) * dDf_dDj;
        pivotsByInf.emplace(infeas, pivotId);
//        std::cout << "Infeasibility error at " << Dj << " = " << infeas << " - " << infeasibility(Dj) << " = "
//        << infeas - infeasibility(Dj) << std::endl;
//        assert(infeas == infeasibility(Dj));
        lastDj = Dj;
        if(pivotId < nonZeroRows.size() * 2) {
            dDf_dDj += fabs(col[nonZeroRows[pivotId / 2]]);
        } else {
            dDf_dDj += 1.0;
        }
    }
    return pivotsByInf;
}

// returns gradient after perturbation of this column by deltaj
double ProposalPivot::colInfeasibilityGradient(double deltaj) {
    double grad = 0.0;
    for(int nzi : nonZeroRows) {
        int rowk = simplex.head[nzi];
        grad += col[nzi] * infeasibilityGradient(simplex.beta[nzi] + col[nzi] * deltaj,
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


void ProposalPivot::setCol(int j) {
    this->j = j;
    col = simplex.tableauCol(j);
    initNonZeroRows();
}

void ProposalPivot::clearCol() {
    for(int i: nonZeroRows) col[i] = 0.0;
    nonZeroRows.clear();
}

// returns true if pmfIndex corresponds to a row that is a structural var and has a unity coefficient
bool ProposalPivot::isActive(int pmfIndex) {
    if(pmfIndex >= 2*nonZeroRows.size()) return true; // bound swaps active (including null pivot)
    // return simplex.isAtUpperBound(j) ^ (pmfIndex%2); // null pivot inactive
    int i = nonZeroRows[pmfIndex / 2];
    if(fabs(fabs(col[i])-1.0) > tol) return false;         // only pivot on unity elements
    int k = simplex.head[i];
    return !simplex.isAuxiliary(k); // don't pivot on auxiliary vars
}

void ProposalPivot::setToPivotIndex(int pivotIndex) {
    leavingVarToUpperBound = pivotIndex % 2;
    if(pivotIndex < 2 * nonZeroRows.size()) {
        i = nonZeroRows[pivotIndex / 2];
        int leavingk = simplex.head[i];
        deltaj = ((leavingVarToUpperBound ? simplex.u[leavingk] : simplex.l[leavingk]) - simplex.beta[i]) / col[i];
        assert(fabs(col[i]) > 1e-7);
        assert(fabs(deltaj) < 10.0);
    } else {
        i = -1; // column does bound swap.
        int k = simplex.head[simplex.nBasic() + j];
        deltaj = (leavingVarToUpperBound?simplex.u[k]:simplex.l[k]) - (simplex.isAtUpperBound(j)?simplex.u[k]:simplex.l[k]);
    }
}

std::vector<double> ProposalPivot::infeasibilityCost() {
    // set objective to out-of-bounds rows
    std::vector<double> infeasibility(simplex.nBasic() + 1);
    infeasibility[0] = 0.0;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        infeasibility[i] = infeasibilityGradient(simplex.beta[i], simplex.l[k], simplex.u[k]);
    }
    return infeasibility;
}


