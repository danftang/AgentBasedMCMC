//
// Created by daniel on 06/07/2021.
//

#ifndef GLPKTEST_POTENTIALENERGYPIVOT_H
#define GLPKTEST_POTENTIALENERGYPIVOT_H

#include <deque>
#include "Phase1Pivot.h"

class PotentialEnergyPivot: public ProposalPivot {
public:
    static constexpr double kappaRow = -7.0;//-6.0;//-3.8;//-1.125; // exponential coefficient for probabilities of choosing row based on change in infeasibility
//    static constexpr double p1 = 0.25; // Given a simplex state that has only one high energy column, p1 gives the probability that the high energy col will be proposed
    double kappaCol = 1.5;//1.5;//1.125;       // exponential coefficient for relative probability of proposing a col based on potential energy set to max(1,log(p1*nNonBasic))
//    double kappaBasis = 3.0; // exponential coefficient for bias towards a particular basis;
//    static constexpr double p0 = 0.01; // relative probability of choosing column with zero reduced cost compared to a high potential col
//    static constexpr double p1 = 0.1; // relative probability of choosing a column with a low potential compared to a high potential col

    // Cached diagnostics for current and proposed states
//    std::vector<double> reducedCost;
//    int infeasibilityCount;
//    double Ep; // current potential energy (minus chosen column's contribution)
    std::vector<double> currentReducedCosts;         // reduced infeasibility cost by column
    std::vector<double> currentPotentialEnergies;    // potential energies by column
    std::vector<double> cdf;                // cumulative distribution function of probabilities of choosing column j

    // These variables contain a cache of costs for the current and proposed states.
    // The first element is for the current value and the last element is for the proposed value.
    // If there's only one element then the current and proposed values are the same.
    std::vector<double> currentInfeasibilityCosts;   // cost by row, depending on infeasibility of each row.
//    std::deque<double>              totalPotentials;
//    std::deque<double>              energyNormalisation;    // sum of exponentials of potential energy

    explicit PotentialEnergyPivot(SimplexMCMC &simplex);

    const PotentialEnergyPivot &nextProposal();

    void initCache();

private:

    bool isInPredPreyPreferredBasis(int j);

    void recalculateCDF();
    void chooseCol();
    void chooseRow();

    void calcAcceptanceContrib();
//    void calcProposedInfeasibilityCosts();
//    void calcProposedReducedCosts();
//    void calcProposedEnergies();

    bool recalculateInfeasibilityCost();
    void recalculateReducedCosts();
    void recalculatePotentials();

//    double potentialEnergy(int j, const std::vector<double> &reducedCost);

    static double potentialEnergy(bool isAtUpperBound, double reducedCost) {
//        if(reducedCost < -0.001) return isAtUpperBound ? 0.0 : 1.0;
//        if(reducedCost > 0.001) return isAtUpperBound ? 1.0 : 0.0;
//        return 0.0;
        return reducedCost * ((isAtUpperBound ? 1.0 : 0.0) - (reducedCost<0.0? 1.0 : 0.0));
    }

    static double infeasibilityCostFn(double b, double lowerBound, double upperBound) {
        if(b < lowerBound - tol) return -1.0;
        if(b > upperBound + tol) return 1.0;
//        if(b < lowerBound + tol) return -0.5;
//        if(b > upperBound - tol) return 0.5;
        return 0.0;
    }

    // debug
//    void checkCurrentCacheIsValid();


//    static int sign(double x);



//    void calcAcceptanceContribCheck();
//    void calcAcceptanceContribCheck();

};


#endif //GLPKTEST_POTENTIALENERGYPIVOT_H
