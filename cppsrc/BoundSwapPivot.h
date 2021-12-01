//
// Created by daniel on 27/11/2021.
//

#ifndef GLPKTEST_BOUNDSWAPPIVOT_H
#define GLPKTEST_BOUNDSWAPPIVOT_H

#include <vector>
#include "glpkpp.h"
#include "ProposalPivot.h"

class BoundSwapPivot: public ProposalPivot {
public:
    static constexpr double kappaRow = -7.5;//-6.5;//-6.0;//-3.8;//-1.125; // exponential coefficient for probabilities of choosing row based on change in infeasibility
    static constexpr double kappaCol = -0.75*kappaRow;//1.5;//1.125;       // exponential coefficient for relative probability of proposing a col based on potential energy set to max(1,log(p1*nNonBasic))

    std::vector<double> currentReducedCosts;         // reduced infeasibility cost by column
    std::vector<double> currentPotentialEnergies;    // potential energies by column
    std::vector<double> cdf;                // cumulative distribution function of probabilities of choosing column j
    std::vector<double> currentInfeasibilityCosts;   // cost by row, depending on infeasibility of each row.

    std::vector<glp::SparseVec> tableauCols;        // The simplex tableau coefficients for this basis

    explicit BoundSwapPivot(SimplexMCMC &simplex);

    const BoundSwapPivot &nextProposal();

    void applyProposal();

    void init();

private:
    void initBasis();
    void calculateTableau();
    bool isInPredPreyPreferredBasis(int k);

    void recalculateCDF();
    void chooseCol();
    void chooseRow();

    void calcAcceptanceContrib();
    bool recalculateInfeasibilityCost();
    void recalculateReducedCosts();
    void recalculatePotentials();


    static double potentialEnergy(bool isAtUpperBound, double reducedCost) {
        if(reducedCost < -0.001) return isAtUpperBound ? 0.0 : 1.0;
        if(reducedCost > 0.001) return isAtUpperBound ? 1.0 : 0.0;
        return 0.0;
    }

    static double infeasibilityCostFn(double b, double lowerBound, double upperBound) {
        if(b < lowerBound - tol) return -1.0;
        if(b > upperBound + tol) return 1.0;
        return 0.0;
    }

    static double infeasibility(double b, double lowerBound, double upperBound) {
        if(b > upperBound) return b - upperBound;
        if(b < lowerBound) return lowerBound - b;
        return 0.0;
    }
};


#endif //GLPKTEST_BOUNDSWAPPIVOT_H
