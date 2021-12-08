//
// Created by daniel on 27/11/2021.
//

#ifndef GLPKTEST_BOUNDSWAPPIVOT_H
#define GLPKTEST_BOUNDSWAPPIVOT_H

#include <vector>
#include "glpkpp.h"
#include "ProposalPivot.h"
#include "MutableCategoricalArray.h"

class BoundSwapPivot {
public:
    static constexpr double kappa = -6.5; // exponential coefficient for probabilities of choosing row based on change in infeasibility

    const int   i;
    int         j;
    bool        leavingVarToUpperBound;
    double      deltaj;                     // the change in the value of this column on pivoting
    SimplexMCMC &simplex;
    double      logAcceptanceContribution;

    std::vector<double> currentDeltaE;              // change of energy on bound swap by column
    std::vector<double> currentE;                   // energy by row
    std::map<int,double> changedDeltaE;             // new values of deltaE by col, where changed by the current proposal
    MutableCategoricalArray colPMF;

    std::vector<glp::SparseVec> tableauCols;        // The simplex tableau coefficients for this basis
    std::vector<glp::SparseVec> tableauRows;
    std::vector<double> feasibleEnergy;             // Energy change in feasible region by kSim

    explicit BoundSwapPivot(SimplexMCMC &simplex);

    BoundSwapPivot &nextProposal();

    void init();

    const glp::SparseVec &tableauCol() const { return tableauCols[j]; }

private:
    void initBasis();
    void calculateTableau();
    bool isInPredPreyPreferredBasis(int k);

    void calcDeltaEChanges();

    double deltaFeasibleEnergy();
    double deltaInfeasibility();

    void chooseCol();
    void chooseBound();

    std::vector<double> calcDeltaE();
    void updateAllCosts();
    void randomiseBounds();

    // debug
    void checkCosts();

//    double infeasibilityCostFn(int i);

//    static double potentialEnergy(bool isAtUpperBound, double reducedCost) {
//        if(reducedCost < -0.001) return isAtUpperBound ? 0.0 : 1.0;
//        if(reducedCost > 0.001) return isAtUpperBound ? 1.0 : 0.0;
//        return 0.0;
//    }

//    static double infeasibilityCostFn(double b, double lowerBound, double upperBound) {
//        if(b < lowerBound - tol) return -1.0;
//        if(b > upperBound + tol) return 1.0;
////        if(b > upperBound - tol) return 0.25; // TODO: test!!!
//        return 0.0;
//    }

    // energy of ith row at beta = b
    double energy(int i, double b);

    static double energy(double b, double feasibleE, double lowerBound, double upperBound) {
        if(b > upperBound) return b - upperBound + feasibleE*(upperBound-lowerBound);
        if(b < lowerBound) return lowerBound - b;
        return (b-lowerBound) * feasibleE;
    }

    static double infeasibility(double b, double lowerBound, double upperBound) {
        if(b > upperBound) return b - upperBound;
        if(b < lowerBound) return lowerBound - b;
        return 0.0;
    }


};


#endif //GLPKTEST_BOUNDSWAPPIVOT_H
