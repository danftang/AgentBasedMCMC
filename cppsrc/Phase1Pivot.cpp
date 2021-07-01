//
// Created by daniel on 28/06/2021.
//

#include "Phase1Pivot.h"
#include "Random.h"
#include <float.h>
#include <cassert>

Phase1Pivot::Phase1Pivot(glp::Simplex &simplex): ProposalPivot(simplex) {
    chooseCol();
    chooseRow();
}

// Choose from among the column(s) with maximum the "potential energy"
// where the "potential energy" of a column j is zero if a perturbation
// of j in the feasible direction would increase the current linearised
// infeasibility, or the absolute value of the gradient if a feasible
// perturbation would decrease the infeasibility.
void Phase1Pivot::chooseCol() {
    setSimplexToInfeasibilityObjective();
    simplex.recalculatePi();

//    // choose random colum with maximum potential
//    std::vector<int> bestCols;
//    double bestGrad = 1.0;
//    for(int j=1; j <= simplex.nNonBasic(); ++j) {
//        double dj = simplex.reducedCost(j);
//        if (dj > tol) {
//            if (simplex.isAtUpperBound(j)) {
//                if (dj > -bestGrad) bestCols.clear();
//                if (dj >= -bestGrad) {
//                    bestGrad = -dj;
//                    bestCols.push_back(j);
//                }
//            }
//        } else if (dj < -tol) {
//            if (!simplex.isAtUpperBound(j)) {
//                if (dj < bestGrad) bestCols.clear();
//                if (dj <= bestGrad) {
//                    bestGrad = dj;
//                    bestCols.push_back(j);
//                }
//            }
//        }
//    }
//    setCol(bestCols[Random::nextInt(0, bestCols.size())]);

    // choose random column with positive potential
    std::vector<int> improvingCols;
    for(int j=1; j <= simplex.nNonBasic(); ++j) {
        if(simplex.isAtUpperBound(j)) {
            if(simplex.reducedCost(j) > tol) {
                improvingCols.push_back(j);
            }
        } else {
            if(simplex.reducedCost(j) < -tol) {
                improvingCols.push_back(j);
            }
        }
    }
    setCol(improvingCols[Random::nextInt(0, improvingCols.size())]);

//    // choose the first column with maximum potential
//    int bestCol = -1;
//    double bestGrad = 1.0;
//    for(int j=1; j <= simplex.nNonBasic(); ++j) {
//        double objj = simplex.reducedObjective(j);
//        if(simplex.isAtUpperBound(j)) {
//            if(objj > tol && objj > -bestGrad) {
//                bestGrad = -objj;
//                bestCol = j;
//            }
//        } else {
//            if(objj < -tol && objj < bestGrad) {
//                bestGrad = objj;
//                bestCol = j;
//            }
//        }
//    }
//    setCol(bestCol);

}


// choose from among the pivots that minimise the infeasibility + delta * the potential energy
// of this column after the pivot (i.e. if there are multiple pivots that minimise infeasibility,
// choose from among the ones that also minimise potential energy of this column).
//
// After a pivot, the gradient of this column is divided by the value of the pivot element,
// so the sign of the gradient changes if the pivot element is negative.
//  So, given a pivot index, p, and assuming dF_dXj != 0, the col is at high potential after the pivot if
//  (p%2) xor (dF_dXj/Mij < 0) = (p%2) xor (dF_dXj < 0) xor (Mij < 0)
// in the case that p < 2*|nonZeroRows| of
// (p%2) xor (dF_dXj < 0)
// otherwise
void Phase1Pivot::chooseRow() {
    std::multimap<double,int> pivots = getPivotsByDeltaJ();
//    assert(simplex.reducedObjective(j) == colInfeasibilityGradient(0.0)); // TEST

    // Choose from among minimum infeasibility
//    static constexpr double deltaPotential = 0.0;
    bool djNegative = simplex.reducedCost(j) < 0.0;
    auto pivotIt = pivots.begin();
    double lastDj = pivotIt->first;
    std::vector<int> bestPivotIndices;
    double infeas = infeasibility(lastDj);
    double minScore = DBL_MAX;
    double dDf_dDj = colInfeasibilityGradient(lastDj - tol);
    for(auto [Dj, pivotIndex] : pivots) {
        infeas += (Dj - lastDj) * dDf_dDj;
        lastDj = Dj;
//        double score = infeas;
//        bool highPotential;
        if(pivotIndex < 2*nonZeroRows.size()) {
            double Mij = col[nonZeroRows[pivotIndex/2]];
//            highPotential = djNegative ^ (pivotIndex % 2) ^ (Mij < 0);
            dDf_dDj += fabs(Mij);
        } else {
//            highPotential = djNegative ^ (pivotIndex % 2);
            dDf_dDj += 1.0;
        }
        if(infeas < minScore) {
            if(isActive(pivotIndex)) {
                bestPivotIndices.clear();
                minScore = infeas;
                bestPivotIndices.push_back(pivotIndex);
            }
//        } else if(infeas == minScore && (!highPotential || Random::nextDouble() > 0.5) && isActive(pivotIndex)) {
        } else if(infeas == minScore && isActive(pivotIndex)) {
            bestPivotIndices.push_back(pivotIndex);
        }
    }
    setToPivotIndex(bestPivotIndices[Random::nextInt(0, bestPivotIndices.size())]);


//    // choose minimum infeasibility with exponentially decreasing prob in the ordering of 'pivots'.
//    auto pivotIt = pivots.begin();
//    double lastDj = pivotIt->first;
//    int minPivotIndex;// = pivotIt->second;
//    double infeas = infeasibility(lastDj);
//    double minScore = DBL_MAX;
//    double dDf_dDj = colInfeasibilityGradient(lastDj - tol);
//    for(auto [Dj, pmfIndex] : pivots) {
//        infeas += (Dj-lastDj)*dDf_dDj;
//        lastDj = Dj;
//        if(pmfIndex < nonZeroRows.size()*2) {
//            dDf_dDj += fabs(col[nonZeroRows[pmfIndex / 2]]);
//        } else {
//            dDf_dDj += 1.0;
//        }
//        if(infeas < minScore) {
//            if(isActive(pmfIndex)) {
//                minScore = infeas;
//                minPivotIndex = pmfIndex;
//            }
//        } else if(infeas == minScore && Random::nextDouble() > 0.5 && isActive(pmfIndex)) {
//            minPivotIndex = pmfIndex;
//        }
//    }
//    setToPivotIndex(minPivotIndex);


    //    while(infeas <= minScore && ++pivotIt != pivots.end()) {
//        auto [Dj, pivotIndex] = *pivotIt;
//        infeas += (Dj-lastDj)*dDf_dDj;
//        assert(infeas == infeasibility(Dj));
//        lastDj = Dj;
//        if(pivotIndex < nonZeroRows.size() * 2) {
//            dDf_dDj += fabs(col[nonZeroRows[pivotIndex / 2]]);
//        } else {
//            dDf_dDj += 1.0;
//        }
////        if((infeas < minScore) && isActive(pivotIndex)) {
//        if(infeas < minScore) {
//            minScore = infeas;
//            minPivotIndex = pivotIndex;
//        }
//    }
//    setToPivotIndex(minPivotIndex);
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
            simplex.c[k] = -1.0;
            ++infeasibilityCount;
        } else if(simplex.b[i] > simplex.u[k] + tol) {
            simplex.c[k] = 1.0;
            ++infeasibilityCount;
        } else {
            simplex.c[k] = 0.0;
        }
    }
    return infeasibilityCount;
}


// returns gradient after perturbation of this column by deltaj
//double Phase1Pivot::colInfeasibilityGradient(double deltaj) {
//    double grad = 0.0;
//    for(int nzi : nonZeroRows) {
//        int rowk = simplex.head[nzi];
//        grad += col[nzi] * infeasibilityGradient(simplex.b[nzi] + col[nzi] * deltaj,
//                                                 simplex.l[rowk],
//                                                 simplex.u[rowk]);
//    }
//    int colk = simplex.head[simplex.nBasic() + j];
//    double xUpperBound = simplex.u[colk];
//    double xLowerBound = simplex.l[colk];
//    double xj = simplex.isAtUpperBound(j)?xUpperBound:xLowerBound;
//    grad += infeasibilityGradient(xj + deltaj, xLowerBound, xUpperBound);
//    return grad;
//}
//
//
//double Phase1Pivot::infeasibilityGradient(double v, double lowerBound, double upperBound) {
//    if(v > upperBound) {
//        return 1.0;
//    } else if(v < lowerBound) {
//        return -1.0;
//    }
//    return 0.0;
//}

// returns true if pmfIndex corresponds to a row that is a structural var and has a unity coefficient
bool Phase1Pivot::isActive(int pmfIndex) { // bound swap active, null pivot inactive
    if(pmfIndex >= 2*nonZeroRows.size()) return simplex.isAtUpperBound(j) ^ (pmfIndex%2);
    int i = nonZeroRows[pmfIndex / 2];
    if(fabs(fabs(col[i])-1.0) > tol) return false;         // only pivot on unity elements
    int k = simplex.head[i];
    return simplex.kSimTokProb[k] > simplex.originalProblem.nConstraints(); // don't pivot on auxiliary vars
}


// returns the infeasibility of the current solution perturbed by this column changing by deltaj
double Phase1Pivot::infeasibility(double deltaj) {
    double dist = 0.0;
    for(int i=1; i <= simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        double v = simplex.b[i] + col[i]*deltaj;
        if(v < simplex.l[k]) {
            dist += simplex.l[k] - v;
        } else if(v > simplex.u[k]) {
            dist += v - simplex.u[k];
        }
    }

    int k = simplex.head[simplex.nBasic() + j];
    double Xj = (simplex.isAtUpperBound(j)?simplex.u[k]:simplex.l[k]) + deltaj;
    if(Xj < simplex.l[k]) {
        dist += simplex.l[k] - Xj;
    } else if(Xj > simplex.u[k]) {
        dist += Xj - simplex.u[k];
    }
    return dist;
}
