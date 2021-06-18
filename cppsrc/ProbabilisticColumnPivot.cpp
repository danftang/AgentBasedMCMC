//
// Created by daniel on 08/06/2021.
//

#include <cmath>
#include <cfloat>
#include "ProbabilisticColumnPivot.h"
#include "Random.h"

void ProbabilisticColumnPivot::chooseCol() {

    // set objective to out-of-bounds rows
    bool isFeasible = true;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        if(simplex.b[i] < simplex.l[k] - tol) {
            simplex.c[i] = -1.0;
            isFeasible = false;
        } else if(simplex.b[i] > simplex.u[k] + tol) {
            simplex.c[i] = 1.0;
            isFeasible = false;
        } else {
            simplex.c[i] = 0.0;
        }
    }
    simplex.recalculatePi();


    if(isFeasible) {
        // null objective so just choose with uniform prob
        j = Random::nextInt(1,simplex.nNonBasic() + 1);
    } else {
        // choose with prob proportional to exp of reduced objective
        std::vector<double> cdf(simplex.nNonBasic() + 1, 0.0);
        std::vector<double> reducedObjective = simplex.reducedObjective();
        double cumulativeP = 0.0;

        for(int j=1; j<= simplex.nNonBasic(); ++j) {
            cumulativeP += reducedObjective[j]<0.0?1.0:0.0;
            cdf[j] = cumulativeP;
        }
        double *it = std::lower_bound(
                cdf.data(),
                cdf.data() + simplex.nNonBasic() + 1,
                Random::nextDouble(0.0, cdf[simplex.nNonBasic()])
                );
        j = it - cdf.data();
    }

}


void ProbabilisticColumnPivot::chooseRow() {

    double feas = feasibility(0.0);
    std::cout << "Feasibility = " << feas << std::endl;

    std::multimap<double, int> transitions; // from delta_j to PMF-index.

    // identify active rows: (rows with non-zero coefficient)
    for(int i=1; i<col.size(); ++i) {
        if(fabs(col[i]) > tol) nonZeroRows.push_back(i);
    }

    // calculate delta_js for each possible pivot
    for(int m=0; m < nonZeroRows.size(); ++m) {
        int i = nonZeroRows[m];
        int k = simplex.head[i];
        double rowLowerBound = simplex.l[k];
        double rowUpperBound = simplex.u[k];
        if(rowLowerBound > -DBL_MAX) {
            double deltaLB = (rowLowerBound - simplex.b[i]) / col[i];
            transitions.emplace(deltaLB, 2 * m);
        }
        if(rowUpperBound < DBL_MAX) {
            double deltaUB = (rowUpperBound - simplex.b[i]) / col[i];
            transitions.emplace(deltaUB, 2 * m + 1);
        }
    }

    // add entry for bounds of this column
    int colk = simplex.head[simplex.nBasic() + j];
    double boundSwapDelta = simplex.isAtUpperBound(j)?(simplex.l[colk]- simplex.u[colk]):(simplex.u[colk]- simplex.l[colk]);
    transitions.emplace(boundSwapDelta, 2 * nonZeroRows.size());
    transitions.emplace(0.0, 2 * nonZeroRows.size() + 1);

    // now populate pivotPMF
    pivotPMF.resize(nonZeroRows.size() * 2 + 1, 0.0);

    double lastDj = transitions.begin()->first;
    double DeltaF = feasibility(lastDj);
    double dDf_dDj = colFeasibilityGradient(lastDj - tol);
    for(auto [Dj, pmfIndex] : transitions) {
        DeltaF += (Dj-lastDj)*dDf_dDj;
        lastDj = Dj;
        if(isActive(pmfIndex)) pivotPMF[pmfIndex] = exp(kappa * DeltaF);
//        std::cout << "delta = " << Dj << " pivotProb = " << pivotPMF[pmfIndex] << " DeltaF = " << DeltaF << " dDf_dDj = " << dDf_dDj << " Feasibility = " << feasibility(Dj) << std::endl;
        if(pmfIndex < nonZeroRows.size()*2) {
            dDf_dDj += fabs(col[nonZeroRows[pmfIndex / 2]]);
        } else {
            dDf_dDj += 1.0;
        }
    }

    // choose row
    int pivotChoice = Random::choose(pivotPMF.begin(), pivotPMF.end());
    if(pivotChoice < 2 * nonZeroRows.size()) {
        i = nonZeroRows[pivotChoice / 2];
        leavingVarToUpperBound = pivotChoice % 2;
        int leavingk = simplex.head[i];
        delta = ((leavingVarToUpperBound?simplex.u[leavingk]:simplex.l[leavingk]) - simplex.b[i])/col[i];
    } else {
        i = -1; // column does bound swap.
        leavingVarToUpperBound = !simplex.isAtUpperBound(j);
        delta = boundSwapDelta;
    }
//    std::cout << "Chose pivot with prob " << pivotPMF[pivotChoice] << " Delta = " << delta << std::endl;
    logTransitionRatio = 0.0;

}


// returns the gradient of the feasibility function with respect to changes in the value of
// this column in either the forward (increasing val) or backward (decreasing val) directions.
// (since the col is on a boundary, the gradient is discontinuous at this point so direction
// must be specified).
//double ProbabilisticColumnPivot::colFeasibilityGradient(bool forward) {
//    double grad = 0.0;
//    for(int nzi : nonZeroRows) {
//        grad += iFeasibilityGradient(nzi, forward != (col[nzi] > 0.0)) * (-col[nzi]);
//    }
//    if(simplex.isAtUpperBound(j)) { // add gradient of this col
//        if(forward) grad += 1.0;
//    } else {
//        if(!forward) grad -= 1.0;
//    }
//    return grad;
//}

// returns gradient after perturbation of this column by deltaj
double ProbabilisticColumnPivot::colFeasibilityGradient(double deltaj) {
    double grad = 0.0;
    for(int nzi : nonZeroRows) {
        int rowk = simplex.head[nzi];
        grad +=  col[nzi] * feasibilityGradient(simplex.b[nzi] + col[nzi]*deltaj,
                                                 simplex.l[rowk],
                                                 simplex.u[rowk]);
    }
    int colk = simplex.head[simplex.nBasic() + j];
    double xUpperBound = simplex.u[colk];
    double xLowerBound = simplex.l[colk];
    double xj = simplex.isAtUpperBound(j)?xUpperBound:xLowerBound;
    grad += feasibilityGradient(xj + deltaj, xLowerBound, xUpperBound);
    return grad;
}


// returns the gradient of the feasibility function for the current value of the i'th row with
// respect to the value of this row.
// If 'forward' is true, returns the gradient in the forward direction, otherwise returns the grad
// in the backward direction (these can be different at boundaries since the gradient
// is discontinuous there).
//double ProbabilisticColumnPivot::iFeasibilityGradient(int i, bool forward) {
//    int k = simplex.head[i];
//    double v = simplex.b[i];
//    if(double upperBound = simplex.u[k]; v > upperBound - forward ? tol : -tol) {
//        return 1.0;
//    } else if(double lowerBound = simplex.l[k]; v < lowerBound - forward ? tol : -tol) {
//        return -1.0;
//    }
//    return 0.0;
//}

double ProbabilisticColumnPivot::feasibilityGradient(double v, double lowerBound, double upperBound) {
    if(v >= upperBound) {
        return 1.0;
    } else if(v < lowerBound) {
        return -1.0;
    }
    return 0.0;
}


//// returns the feasibility of the current solution in simplex
//double ProbabilisticColumnPivot::feasibility(int j) {
//    double dist = 0.0;
//    for(int i=1; i < simplex.nBasic(); ++i) {
//        int k = simplex.head[i];
//        if(simplex.b[i] < simplex.l[k]) {
//            dist += simplex.l[k] - simplex.b[i];
//        } else if(simplex.b[i] > simplex.u[k]) {
//            dist += simplex.b[i] - simplex.u[k];
//        }
//    }
//    return dist;
//}

// returns the feasibility of the current solution perturbed by this column changing by deltaj
double ProbabilisticColumnPivot::feasibility(double deltaj) {
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
    if(simplex.isAtUpperBound(j)) {
        if(deltaj > 0) dist += deltaj;
    } else {
        if(deltaj < 0) dist -= deltaj;
    }
    return dist;
}


// returns true if pmfIndex corresponds to a row that is a structural var and has a unity coefficient
bool ProbabilisticColumnPivot::isActive(int pmfIndex) {
    if(pmfIndex == 2 * nonZeroRows.size()) return true;         // bound swap always active
    if(pmfIndex == 2 * nonZeroRows.size() + 1) return false;    // null pivot never active
    int i = nonZeroRows[pmfIndex / 2];
    if(fabs(fabs(col[i])-1.0) > tol) return false;         // only pivot on unity elements
    int k = simplex.head[i];
    return simplex.kSimTokProb[k] > simplex.originalProblem.nConstraints(); // don't pivot on auxiliary vars
}
