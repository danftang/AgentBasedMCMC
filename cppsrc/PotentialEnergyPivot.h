//
// Created by daniel on 06/07/2021.
//

#ifndef GLPKTEST_POTENTIALENERGYPIVOT_H
#define GLPKTEST_POTENTIALENERGYPIVOT_H


#include "glpkpp.h"
#include "Phase1Pivot.h"

class PotentialEnergyPivot: public Phase1Pivot {
public:
    static constexpr double kappaRow = -1.0;// 1.0;//-5.0; // exponential coefficient for probabilities of choosing row based on change in infeasibility
    static constexpr double p1 = 0.2; // Given a simplex state that has only one high energy column, p1 gives the probability that the high energy col will be proposed
    double kappaCol;        // exponential coefficient for relative probability of proposing a col based on potential energy set to max(1,log(p1*nNonBasic))
//    static constexpr double p0 = 0.01; // relative probability of choosing column with zero reduced cost compared to a high potential col
//    static constexpr double p1 = 0.1; // relative probability of choosing a column with a low potential compared to a high potential col

    std::vector<double> reducedCost;
    int infeasibilityCount;
    double Ep; // current potential energy (minus chosen column's contribution)

    explicit PotentialEnergyPivot(glp::Simplex &simplex);

    void chooseCol();
    void chooseRow();

    void calcAcceptanceContrib();
    double potentialEnergy(int j, const std::vector<double> &reducedCost) {
//        return sign(reducedCost[j]) * (simplex.isAtUpperBound(j) ? 1.0 : 0.0);
//        return reducedCost[j] * (simplex.isAtUpperBound(j) ? 1.0 : 0.0);
        return reducedCost[j] * ((simplex.isAtUpperBound(j) ? 1.0 : 0.0) - (reducedCost[j]<0.0? 1.0 : 0.0));
    }

    static int sign(double x);

    void calcAcceptanceContrib2();
};


#endif //GLPKTEST_POTENTIALENERGYPIVOT_H
