//
// Created by daniel on 06/07/2021.
//

#include <cfloat>
#include "PotentialEnergyPivot.h"
#include "Random.h"

PotentialEnergyPivot::PotentialEnergyPivot(glp::Simplex &simplex):
Phase1Pivot(simplex,0,0) {
    chooseCol();
    chooseRow();
    calcAcceptanceContrib();
}


// Chooses columns based on the exponential of their "Potential energy" which is defined as d_j*(2X_j-1)
// where d_j is the reduced cost of the column and X_j is its value, and a column is chosen with
// a probaqbility equal to e^kEj/sum_l e^kEl, where Ej is the potential energy of column j.
//
// The porobability of a state is penalised by A * sum_l e^-kEl * e^-kE / n where E is the total potential energy
// summed over all non-basic columns and n is the number of non-basics. E is observed to be non-negative
// (TODO: should we limit penalty to be always less than 1?).
//
// As a consequence, the contribution of the column choice to the acceptance probability is
// e^-k(DE - DEj) where DE and DEj are the change in total and column energy respectively, during the pivot.
void PotentialEnergyPivot::chooseCol() {
    // set objective to out-of-bounds rows
    infeasibilityCount = setSimplexToInfeasibilityObjective();

    if(infeasibilityCount == 0) {
        // null objective so just chooseFromPMF with uniform prob
        setCol(Random::nextInt(1,simplex.nNonBasic() + 1));
    } else {
        // chooseFromPMF with prob proportional to exponential of potential energy
        simplex.recalculatePi();
        reducedCost = simplex.reducedCost();
        std::vector<double> cdf(simplex.nNonBasic() + 1, 0.0);
        double cumulativeP = 0.0;
        for(int j=1; j<= simplex.nNonBasic(); ++j) {
            double colPotential = reducedCost[j] * (simplex.isAtUpperBound(j) ? 1.0 : -1.0);
            if(colPotential > tol) colPotential = 1.0; else colPotential = 0.0; // TODO: TEST
            cumulativeP += exp(kappaCol*colPotential);
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



void PotentialEnergyPivot::chooseRow() {

    ///////////////////////// debug only
    // double infeas = infeasibility(0.0);
    //   if(infeas > 0.0) std::cout << "Infeasibility = " << infeas << std::endl;
    //////////////

    std::multimap<double, int> transitions = getPivotsByDeltaJ(); // from delta_j to PMF-index.

    // now populate pivotPMF by going from lowest to highest delta_j in order
    std::vector<double> pivotPMF(nonZeroRows.size() * 2 + 2, DBL_MAX); // index is (2*activeRowIndex + toUpperBound), value is probability mass
    double lastDj = transitions.begin()->first;
    double feas = 0.0; //infeasibility(lastDj) - infeasibility(0.0);
    double feasMin = DBL_MAX;
    double dDf_dDj = colInfeasibilityGradient(lastDj - tol);
    for(auto [Dj, pmfIndex] : transitions) {
        feas += (Dj - lastDj) * dDf_dDj;
        lastDj = Dj;
        if(isActive(pmfIndex)) {
            pivotPMF[pmfIndex] = feas;
            if(feas < feasMin) feasMin = feas;
        }
//        std::cout << "deltaj = " << Dj << " pivotProb = " << pivotPMF[pmfIndex] << " feas = " << feas << " dDf_dDj = " << dDf_dDj << " Feasibility = " << infeasibility(Dj) << std::endl;
        if(pmfIndex < nonZeroRows.size()*2) {
            dDf_dDj += fabs(col[nonZeroRows[pmfIndex / 2]]);
        } else {
            dDf_dDj += 1.0;
        }
    }
    assert(feasMin != DBL_MAX);

    // rescale and take exponential
    for(int i=0; i<pivotPMF.size(); ++i) {
        pivotPMF[i] = exp(kappaRow*(pivotPMF[i] - feasMin));
    }
    pivotPMF[2*nonZeroRows.size() + simplex.isAtUpperBound(j)] *= 0.1; // preference against null pivot if other possibilities exist

    // chooseFromPMF row

    setToPivotIndex(Random::chooseFromPMF(pivotPMF.begin(), pivotPMF.end()));

    //    std::cout << "Chose pivot with prob " << pivotPMF[pivotChoice] << " Delta = " << deltaj << std::endl;
//  std::cout << "Proposing pivot from " << !sourceObjectiveIsZero << " to " << destinationObjective << std::endl;
}


// The acceptance contribution is given by
// A = e^-kCol(DE - DEj)
// where DE is the change in potential energy and DEj is the change in potential energy of the pivot column.
// We make use of the fact that the change in the reduced cost of column q!=j is given by
// Ddq = Tiq (dj / Tij)
// where Ti is the pivot row in the tableau.
// And, since the values of the columns q!=j do not change
// DE-DEj = (dj / -Tij) sum_{q!=j} Tiq (2Xq - 1.0)
//
// In the case of a bound swap, the reduced cost does not change so DE-DEj = 0.
void PotentialEnergyPivot::calcAcceptanceContrib() {
//    logAcceptanceContribution = 0.0;
    if(i <= 0 || infeasibilityCount == 0 || reducedCost[j] == 0.0) {
        logAcceptanceContribution = 0.0;
    } else {
        if(deltaj != 0.0) {
            // TODO: Calculate contribution non-degenerate pivots
            logAcceptanceContribution = 0.0;
        } else {
            std::vector<double> Ti = simplex.tableauRow(i);
            double DEnergy = 0.0;
            double dj_Tij = reducedCost[j] / Ti[j];
            for (int q = 1; q <= simplex.nNonBasic(); ++q) {
                if (q != j) {
                    double newdq = reducedCost[q] - dj_Tij * Ti[q];
                    double DEq = ((newdq > tol) - (reducedCost[q] > tol)) * (simplex.isAtUpperBound(q) ? 1.0 : -1.0);
                    DEnergy += DEq;
                }
            }
//        DEnergy *= -reducedCost[j] / Ti[j];
//            std::cout << "DEnergy = " << DEnergy << std::endl;
            logAcceptanceContribution = -kappaCol * DEnergy;
        }
    }
}
