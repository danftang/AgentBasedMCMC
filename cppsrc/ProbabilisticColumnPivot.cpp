//
// Created by daniel on 08/06/2021.
//

#include <cmath>
#include <cfloat>
#include "ProbabilisticColumnPivot.h"
#include "Random.h"



ProbabilisticColumnPivot::ProbabilisticColumnPivot(glp::Simplex &simplex) : simplex(simplex) {
    chooseCol();
    col = simplex.tableauCol(j);
    chooseRow();
}


// Chooses columns based on whether the coefficient of the reduced infeasibility objective
// is zero or not. Non-zero coefficients are alpha times more likely to be chosen than zero coefficients.
//
// The probability of choosing a non-zero coeff (if there is one) is alpha/((alpha-1)*nnz + n), where nnz is the number of
// non-zero coefficients in the reduced objective and n is the total number of non-basics.
// The probability of choosing a zero coeff is 1/((alpha-1)*nnz + n).
// If we also penalise the probability of infeasible states by a factor of A*((alpha-1)*nnz + n)*alpha^-n_inf
// where n_inf is the number of infeasible columns
// then the overall contribution to the acceptance prob is either A*alpha^(1-n_inf) or A*alpha^-n_inf depending on whether
// the column has non-zero and zero coeff respectively.
// So the backwards/forwards ratio is
//  * alpha^(-Dn_Inf) for non-zero -> non-zero and zero -> zero transitions
//  * alpha^(-(Dn_inf-1)) for zero -> non-zero transitions
//  * alpha^(-(Dn_inf+1)) for non-zero -> zero transitions
// But in the case of a non-zero -> zero transition Dn_inf is likely to be negative, so we are likely to get acceptance
// and in the case of zero-> non-zero Dn_inf is likely to be positive so, again, we get acceptance.
// Note, we only need to calculate n_inf and d_j for the source and destination, which we would have to calculate anyway
// for the row calculation.
// Note also that the destination d_j is zero if and only if the current column has d_j of zero under the destination B
// so there's no need to re-calculate the column for the destination.
void ProbabilisticColumnPivot::chooseCol() {
    // set objective to out-of-bounds rows
    sourceInfeasibilityCount = 0;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        if(simplex.b[i] < simplex.l[k] - tol) {
            simplex.c[i] = -1.0;
            ++sourceInfeasibilityCount;
        } else if(simplex.b[i] > simplex.u[k] + tol) {
            simplex.c[i] = 1.0;
            ++sourceInfeasibilityCount;
        } else {
            simplex.c[i] = 0.0;
        }
    }


    if(sourceInfeasibilityCount == 0) {
        // null objective so just choose with uniform prob
        sourceObjectiveIsZero = true;
        j = Random::nextInt(1,simplex.nNonBasic() + 1);
    } else {
        // choose with prob proportional to alpha if reduced objective is non-zero or 1 if zero
        simplex.recalculatePi();
        std::vector<double> cdf(simplex.nNonBasic() + 1, 0.0);
        std::vector<double> reducedObjective = simplex.reducedObjective();
        double cumulativeP = 0.0;

        for(int j=1; j<= simplex.nNonBasic(); ++j) {
            cumulativeP += fabs(reducedObjective[j])>tol?alpha:1.0;
            cdf[j] = cumulativeP;
        }
        double *it = std::lower_bound(
                cdf.data(),
                cdf.data() + simplex.nNonBasic() + 1,
                Random::nextDouble(0.0, cdf[simplex.nNonBasic()])
                );
        sourceObjectiveIsZero = (*it == 1.0); // record this for calculation of acceptance contribution later.
        j = it - cdf.data();
    }

}


void ProbabilisticColumnPivot::chooseRow() {

    ///////////////////////// debug only
    double feas = infeasibility(0.0);
    std::cout << "Feasibility = " << feas << std::endl;
    //////////////

    std::multimap<double, int> transitions; // from delta_j to PMF-index.

    // identify active rows: (rows with non-zero coefficient)
    for(int i=1; i<col.size(); ++i) {
        if(fabs(col[i]) > tol) nonZeroRows.push_back(i);
    }

    // calculate delta_js for each possible pivot and place in "transitions" to sort
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

    // add transitions entries for bounds of this column
    int colk = simplex.head[simplex.nBasic() + j];
    double boundSwapDelta = simplex.isAtUpperBound(j)?(simplex.l[colk]- simplex.u[colk]):(simplex.u[colk]- simplex.l[colk]);
    transitions.emplace(boundSwapDelta, 2 * nonZeroRows.size());
    transitions.emplace(0.0, 2 * nonZeroRows.size() + 1);

    // now populate pivotPMF by going from lowest to highest delta_j in order
    pivotPMF.resize(nonZeroRows.size() * 2 + 1, 0.0);
    double lastDj = transitions.begin()->first;
    double DeltaF = infeasibility(lastDj);
    double dDf_dDj = colFeasibilityGradient(lastDj - tol);
    for(auto [Dj, pmfIndex] : transitions) {
        DeltaF += (Dj-lastDj)*dDf_dDj;
        lastDj = Dj;
        if(isActive(pmfIndex)) pivotPMF[pmfIndex] = exp(kappa * DeltaF);
//        std::cout << "deltaj = " << Dj << " pivotProb = " << pivotPMF[pmfIndex] << " DeltaF = " << DeltaF << " dDf_dDj = " << dDf_dDj << " Feasibility = " << infeasibility(Dj) << std::endl;
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
        deltaj = ((leavingVarToUpperBound ? simplex.u[leavingk] : simplex.l[leavingk]) - simplex.b[i]) / col[i];
    } else {
        i = -1; // column does bound swap.
        leavingVarToUpperBound = !simplex.isAtUpperBound(j);
        deltaj = boundSwapDelta;
    }
//    std::cout << "Chose pivot with prob " << pivotPMF[pivotChoice] << " Delta = " << deltaj << std::endl;

    // calculate acceptance contribution
    int destinationInfeasibilityCount = 0;
    double destinationObjective = 0.0;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        double bi = simplex.b[i] + deltaj * col[i];
        if(bi < simplex.l[k] - tol) {
            ++destinationInfeasibilityCount;
            destinationObjective -= col[i];
        } else if(bi > simplex.u[k] + tol) {
            ++destinationInfeasibilityCount;
            destinationObjective += col[i];
        }
    }
    int DInfeasibilityCount = destinationInfeasibilityCount - sourceInfeasibilityCount;
    if(sourceInfeasibilityCount == 0) {
        if(destinationInfeasibilityCount != 0) DInfeasibilityCount -= 1;
    } else if(destinationInfeasibilityCount == 0) {
        if(destinationInfeasibilityCount == 0) DInfeasibilityCount += 1;
    }
    logAcceptanceContribution = -DInfeasibilityCount * log(alpha);
}


// returns the gradient of the infeasibility function with respect to changes in the value of
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


// returns the gradient of the infeasibility function for the current value of the i'th row with
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

double ProbabilisticColumnPivot::infeasibilityGradient(double v, double lowerBound, double upperBound) {
    if(v >= upperBound) {
        return 1.0;
    } else if(v < lowerBound) {
        return -1.0;
    }
    return 0.0;
}


//// returns the infeasibility of the current solution in simplex
//double ProbabilisticColumnPivot::infeasibility(int j) {
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

// returns the infeasibility of the current solution perturbed by this column changing by deltaj
double ProbabilisticColumnPivot::infeasibility(double deltaj) {
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


int ProbabilisticColumnPivot::infeasibilityCount(double deltaj) {
    int count = 0;
    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        if(double bi = simplex.b[i] + deltaj*col[i]; (bi < simplex.l[k] - tol) || (bi > simplex.u[k] + tol))
            ++count;
    }
    return count;
}
