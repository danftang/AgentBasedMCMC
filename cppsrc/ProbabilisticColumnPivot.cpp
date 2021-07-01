//
// Created by daniel on 08/06/2021.
//

#include <cmath>
#include <cfloat>
#include "ProbabilisticColumnPivot.h"
#include "Random.h"



ProbabilisticColumnPivot::ProbabilisticColumnPivot(glp::Simplex &simplex) : Phase1Pivot(simplex,0,0) {
    chooseCol();
    chooseRow();
    calcAcceptanceContrib();
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
    int infeasibilityCount = setSimplexToInfeasibilityObjective();

    if(infeasibilityCount == 0) {
        // null objective so just chooseFromPMF with uniform prob
        setCol(Random::nextInt(1,simplex.nNonBasic() + 1));
    } else {
        // chooseFromPMF with prob proportional to alpha if reduced objective is non-zero or 1 if zero
        simplex.recalculatePi();
        std::vector<double> cdf(simplex.nNonBasic() + 1, 0.0);
        std::vector<double> reducedCost = simplex.reducedCost();
        double cumulativeP = 0.0;

        for(int j=1; j<= simplex.nNonBasic(); ++j) {
            cumulativeP += fabs(reducedCost[j]) > tol ? alpha : 1.0;
            cdf[j] = cumulativeP;
        }
        double *it = std::lower_bound(
                cdf.data(),
                cdf.data() + simplex.nNonBasic() + 1,
                Random::nextDouble(0.0, cdf[simplex.nNonBasic()])
                );
        setCol(it - cdf.data());
    }

}



void ProbabilisticColumnPivot::chooseRow() {

    ///////////////////////// debug only
    double infeas = infeasibility(0.0);
    if(infeas > 0.0) std::cout << "Infeasibility = " << infeas << std::endl;
    //////////////

    std::multimap<double, int> transitions = getPivotsByDeltaJ(); // from delta_j to PMF-index.

    // now populate pivotPMF by going from lowest to highest delta_j in order
    std::vector<double> pivotPMF(nonZeroRows.size() * 2 + 2, 0.0); // index is (2*activeRowIndex + toUpperBound), value is probability mass
    double lastDj = transitions.begin()->first;
    double DeltaF = infeasibility(lastDj); // TODO: do we really need to calculate this?
    double dDf_dDj = colInfeasibilityGradient(lastDj - tol);
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

    // chooseFromPMF row
    setToPivotIndex(Random::chooseFromPMF(pivotPMF.begin(), pivotPMF.end()));
//    std::cout << "Chose pivot with prob " << pivotPMF[pivotChoice] << " Delta = " << deltaj << std::endl;

    // test
//    std::cout << "Proposing pivot from " << !sourceObjectiveIsZero << " to " << destinationObjective << std::endl;

}


void ProbabilisticColumnPivot::calcAcceptanceContrib() {
    int destinationInfeasibilityCount = 0;
    int sourceInfeasibilityCount = 0;

    for(int i=1; i<=simplex.nBasic(); ++i) {
        int k = simplex.head[i];
        double bis = simplex.b[i];
        double bid = bis + deltaj * col[i];
        double lowerBound = simplex.l[k] - tol;
        double upperBound = simplex.u[k] + tol;
        if(bis < lowerBound) {
            ++sourceInfeasibilityCount;
        } else if(bis > upperBound) {
            ++sourceInfeasibilityCount;
        }
        if(bid < lowerBound) {
            ++destinationInfeasibilityCount;
        } else if(bid > upperBound) {
            ++destinationInfeasibilityCount;
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



// returns the infeasibility of the current solution perturbed by this column changing by deltaj
//double ProbabilisticColumnPivot::infeasibility(double deltaj) {
//    double dist = 0.0;
//    for(int i=1; i <= simplex.nBasic(); ++i) {
//        int k = simplex.head[i];
//        double v = simplex.b[i] + col[i]*deltaj;
//        if(v < simplex.l[k]) {
//            dist += simplex.l[k] - v;
//        } else if(v > simplex.u[k]) {
//            dist += v - simplex.u[k];
//        }
//    }
//    if(simplex.isAtUpperBound(j)) {
//        if(deltaj > 0) dist += deltaj;
//    } else {
//        if(deltaj < 0) dist -= deltaj;
//    }
//    return dist;
//}


// returns true if pmfIndex corresponds to a row that is a structural var and has a unity coefficient
//bool ProbabilisticColumnPivot::isActive(int pmfIndex) {
//    if(pmfIndex == 2 * nonZeroRows.size()) return true;         // bound swap always active
//    if(pmfIndex == 2 * nonZeroRows.size() + 1) return false;    // null pivot never active
//    int i = nonZeroRows[pmfIndex / 2];
//    if(fabs(fabs(col[i])-1.0) > tol) return false;         // only pivot on unity elements
//    int k = simplex.head[i];
//    return simplex.kSimTokProb[k] > simplex.originalProblem.nConstraints(); // don't pivot on auxiliary vars
//}


//int ProbabilisticColumnPivot::infeasibilityCount(double deltaj) {
//    int count = 0;
//    for(int i=1; i<=simplex.nBasic(); ++i) {
//        int k = simplex.head[i];
//        if(double bi = simplex.b[i] + deltaj*col[i]; (bi < simplex.l[k] - tol) || (bi > simplex.u[k] + tol))
//            ++count;
//    }
//    return count;
//}
