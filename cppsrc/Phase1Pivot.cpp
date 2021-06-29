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


void Phase1Pivot::chooseCol() {
    setSimplexToInfeasibilityObjective();
    simplex.recalculatePi();

    std::vector<int> improvingCols;
    for(int j=1; j <= simplex.nNonBasic(); ++j) {
        if(simplex.isAtUpperBound(j)) {
            if(simplex.reducedObjective(j) > tol) improvingCols.push_back(j);
        } else {
            if(simplex.reducedObjective(j) < -tol) improvingCols.push_back(j);
        }
    }
    setCol(improvingCols[Random::nextInt(0, improvingCols.size())]);

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


void Phase1Pivot::chooseRow() {
    std::multimap<double,int> pivots = getPivotsByDeltaJ();

//    assert(simplex.reducedObjective(j) == colInfeasibilityGradient(0.0)); // TEST

    auto pivotIt = pivots.begin();
    double lastDj = pivotIt->first;
    int minPivotIndex = pivotIt->second;
    double DeltaF = infeasibility(lastDj);
    double minDeltaF = DeltaF;
    double dDf_dDj = colInfeasibilityGradient(lastDj - tol);
    for(auto [Dj, pmfIndex] : pivots) {
        DeltaF += (Dj-lastDj)*dDf_dDj;
//        std::cout << "Infeasibility error at " << Dj << " = " << DeltaF << " - " << infeasibility(Dj) << " = "
//        << DeltaF - infeasibility(Dj) << std::endl;
        assert(DeltaF == infeasibility(Dj));
        lastDj = Dj;
        if(pmfIndex < nonZeroRows.size()*2) {
            dDf_dDj += fabs(col[nonZeroRows[pmfIndex / 2]]);
        } else {
            dDf_dDj += 1.0;
        }
        if(DeltaF < minDeltaF) {
            if(isActive(pmfIndex)) {
                minDeltaF = DeltaF;
                minPivotIndex = pmfIndex;
            }
        } else if(DeltaF == minDeltaF && Random::nextDouble() > 0.5 && isActive(pmfIndex)) {
            minPivotIndex = pmfIndex;
        }
    }

//    while(DeltaF <= minDeltaF && ++pivotIt != pivots.end()) {
//        auto [Dj, pivotIndex] = *pivotIt;
//        DeltaF += (Dj-lastDj)*dDf_dDj;
//        assert(DeltaF == infeasibility(Dj));
//        lastDj = Dj;
//        if(pivotIndex < nonZeroRows.size() * 2) {
//            dDf_dDj += fabs(col[nonZeroRows[pivotIndex / 2]]);
//        } else {
//            dDf_dDj += 1.0;
//        }
////        if((DeltaF < minDeltaF) && isActive(pivotIndex)) {
//        if(DeltaF < minDeltaF) {
//            minDeltaF = DeltaF;
//            minPivotIndex = pivotIndex;
//        }
//    }
    setToPivotIndex(minPivotIndex);
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
    if(v > upperBound) {
        return 1.0;
    } else if(v < lowerBound) {
        return -1.0;
    }
    return 0.0;
}

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
