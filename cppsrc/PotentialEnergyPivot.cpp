//
// Created by daniel on 06/07/2021.
//

#include <cfloat>
#include "PotentialEnergyPivot.h"
#include "Random.h"

PotentialEnergyPivot::PotentialEnergyPivot(glp::Simplex &simplex): Phase1Pivot(simplex,0,0) {
    chooseCol();
    chooseRow();
//    calcAcceptanceContrib();
}


// Chooses columns based on the sign of their "Potential energy" which is defined as d_j*(2X_j-1)
// where d_j is the reduced cost of the column and X_j is its value.
// If the potential energy is positive, choose with a probability porportional to 1
// if negative, choose wth a prob proportional to p1
// if zero, choose with a prob of p0
void PotentialEnergyPivot::chooseCol() {
    // set objective to out-of-bounds rows
    int infeasibilityCount = setSimplexToInfeasibilityObjective();

    if(infeasibilityCount == 0) {
        // null objective so just chooseFromPMF with uniform prob
        setCol(Random::nextInt(1,simplex.nNonBasic() + 1));
    } else {
        // chooseFromPMF with prob depending on sign of potential energy
        simplex.recalculatePi();
        std::vector<double> cdf(simplex.nNonBasic() + 1, 0.0);
        std::vector<double> reducedCost = simplex.reducedCost();
        double cumulativeP = 0.0;

        for(int j=1; j<= simplex.nNonBasic(); ++j) {
            double potential = reducedCost[j] * (simplex.isAtUpperBound(j)?1.0:-1.0);
            if(potential > tol) {
                cumulativeP += 1.0;
            } else if(potential < -tol) {
                cumulativeP += p1;
            } else {
                cumulativeP += p0;
            }
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
    double infeas = infeasibility(0.0);
    //   if(infeas > 0.0) std::cout << "Infeasibility = " << infeas << std::endl;
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
        if(isActive(pmfIndex)) pivotPMF[pmfIndex] = exp(kappa * DeltaF) + DBL_EPSILON;
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


