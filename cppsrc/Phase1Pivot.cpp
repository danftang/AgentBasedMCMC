//
// Created by daniel on 28/06/2021.
//

#include "Phase1Pivot.h"
#include "Random.h"
#include <float.h>
#include <cassert>

Phase1Pivot::Phase1Pivot(glp::Simplex &simplex): ProposalPivot(simplex), infeasibilityGradient(simplex.nBasic()+1) {
    initInfeasibilityGradient();
    simplex.btran(infeasibilityGradient); // turn objective into pi.
    chooseCol();
    chooseRow();
}

// Choose from among the column(s) with maximum the "potential energy"
// where the "potential energy" of a column j is zero if a perturbation
// of j in the feasible direction would increase the current linearised
// infeasibility, or the absolute value of the gradient if a feasible
// perturbation would decrease the infeasibility.
void Phase1Pivot::chooseCol() {
//    initInfeasibilityGradient();
//    simplex.recalculatePi();

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
//    int nnz = 0;
//    double potential = 0.0;
    double sumOfReducedCost = 0.0;
    std::vector<double> reducedCost = simplex.piTimesMinusN(infeasibilityGradient);
    for(int j=1; j <= simplex.nNonBasic(); ++j) {
        double dj = reducedCost[j];
//        potential += dj * (simplex.isAtUpperBound(j)?0.5:-0.5); // TODO: for Debug only
//        potential += dj * (simplex.isAtUpperBound(j)?1.0:0.0); // TODO: for Debug only
        sumOfReducedCost += dj;
        if(dj > tol) {
//            ++nnz;
            if(simplex.isAtUpperBound(j)) improvingCols.push_back(j);
        } else if(dj < -tol) {
//            ++nnz;
            if(!simplex.isAtUpperBound(j)) improvingCols.push_back(j);
        }
    }
    setCol(improvingCols[Random::nextInt(0, improvingCols.size())]);
//    std::cout << "potential ratio = " << improvingCols.size() << " / " << nnz << " = " << improvingCols.size()*improvingCols.size()*1.0/nnz << std::endl;
//    std::cout << improvingCols.size()*improvingCols.size()*1.0/nnz << std::endl; // seems to be a good monotonically(ish) decreasing value
//    std::cout << improvingCols.size()*infeasibility()*1.0/nnz << std::endl; // seems to be a very good monotonically decreasing value
//    std::cout << improvingCols.size()-(nnz-improvingCols.size()) << std::endl; // seems to be a very good monotonically decreasing value not using infeasibility
//    std::cout << improvingCols.size() << " " << nnz << " " << infeasibility() << " " << infeasibilityCount() << std::endl; // seems to be a good monotonically(ish) decreasing value

//    std::cout << nnz << std::endl;
//    std::cout << infeasibility() - 0.01*sumOfReducedCost << " " << potential << std::endl;
//    std::cout << 6*infeasibility() - sumOfReducedCost << std::endl;

//    std::cout << potential - simplex.reducedCost(j) * (simplex.isAtUpperBound(j)?0.5:-0.5) << std::endl;

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

    // Choose from among minimum infeasibility with uniform prob
//    static constexpr double deltaPotential = 0.0;
//    int reducedCostj = simplex.reducedCost(j);
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
            double Tij = col[nonZeroRows[pivotIndex/2]];
//            highPotential = djNegative ^ (pivotIndex % 2) ^ (Mij < 0);
            dDf_dDj += fabs(Tij);
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
//        } else if(infeas == minScore && (highPotential || Random::nextDouble() < 0.1) && isActive(pivotIndex)) {
        } else if(infeas == minScore && isActive(pivotIndex)) {
            bestPivotIndices.push_back(pivotIndex);
        }
    }
    assert(bestPivotIndices.size() != 0);
    setToPivotIndex(bestPivotIndices[Random::nextInt(0, bestPivotIndices.size())]);

//    if(deltaj == 0.0) {
//        std::cout << "Change in potential = " << (reducedCostj/col[i])*(leavingVarToUpperBound?0.5:-0.5) << " - " << reducedCostj*(simplex.isAtUpperBound(j)?0.5:-0.5)
//        << " = " << (reducedCostj/col[i])*(leavingVarToUpperBound?0.5:-0.5) - reducedCostj*(simplex.isAtUpperBound(j)?0.5:-0.5) << std::endl;
//    }
//    if(deltaj==0.0 && i>0 && !((reducedCostj<0) ^ leavingVarToUpperBound ^ (col[i]<0))) std::cout << "chosen low potential leaving var" << std::endl;
//    if(deltaj==0.0 && i<=0 && !((reducedCostj<0) ^ leavingVarToUpperBound)) std::cout << "chosen low potential bound swap" << std::endl;

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
        deltaj = ((leavingVarToUpperBound ? simplex.u[leavingk] : simplex.l[leavingk]) - simplex.beta[i]) / col[i];
        assert(fabs(col[i]) > 1e-7);
        assert(fabs(deltaj) < 10.0);
    } else {
        i = -1; // column does bound swap.
        int k = simplex.head[simplex.nBasic() + j];
        deltaj = (leavingVarToUpperBound?simplex.u[k]:simplex.l[k]) - (simplex.isAtUpperBound(j)?simplex.u[k]:simplex.l[k]);
    }
}


int Phase1Pivot::initInfeasibilityGradient() {
    // set objective to out-of-bounds rows

    int infeasibilityCount = 0;
    infeasibilityGradient[0] = 0.0;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        if(simplex.beta[i] < simplex.l[k] - tol) {
            infeasibilityGradient[i] = -1.0;
            ++infeasibilityCount;
        } else if(simplex.beta[i] > simplex.u[k] + tol) {
            infeasibilityGradient[i] = 1.0;
            ++infeasibilityCount;
        } else {
            infeasibilityGradient[i] = 0.0;
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
bool Phase1Pivot::isActive(int pmfIndex) {
    if(pmfIndex >= 2*nonZeroRows.size()) return true; // bound swaps active (including null pivot)
        // return simplex.isAtUpperBound(j) ^ (pmfIndex%2); // null pivot inactive
    int i = nonZeroRows[pmfIndex / 2];
    if(fabs(fabs(col[i])-1.0) > tol) return false;         // only pivot on unity elements
    int k = simplex.head[i];
    return !simplex.isAuxiliary(k); // don't pivot on auxiliary vars
}


// returns the infeasibility of the current exactEndState perturbed by this column changing by deltaj
double Phase1Pivot::infeasibility(double deltaj) {
    double dist = 0.0;
    for(int i=1; i <= simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        double v = simplex.beta[i] + col[i]*deltaj;
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

double Phase1Pivot::infeasibility() {
    double dist = 0.0;
    for(int i=1; i <= simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        double v = simplex.beta[i];
        if(v < simplex.l[k]) {
            dist += simplex.l[k] - v;
        } else if(v > simplex.u[k]) {
            dist += v - simplex.u[k];
        }
    }
    return dist;
}

int Phase1Pivot::infeasibilityCount() {
    int infeasibility = 0;
    for(int i=1; i<=simplex.nBasic();++i) {
        int k = simplex.head[i];
        if(simplex.beta[i] < simplex.l[k] - tol || simplex.beta[i] > simplex.u[k] + tol) ++infeasibility;
    }
    return infeasibility;
}
