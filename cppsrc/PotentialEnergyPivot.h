//
// Created by daniel on 06/07/2021.
//

#ifndef GLPKTEST_POTENTIALENERGYPIVOT_H
#define GLPKTEST_POTENTIALENERGYPIVOT_H


#include "glpkpp.h"
#include "Phase1Pivot.h"

class PotentialEnergyPivot: public Phase1Pivot {
public:
    static constexpr double kappa = -5.0; // exponential coefficient for probabilities of choosing row based on change in infeasibility
    static constexpr double p0 = 0.01; // relative probability of choosing column with zero reduced cost compared to a high potential col
    static constexpr double p1 = 0.1; // relative probability of choosing a column with a low potential compared to a high potential col

    explicit PotentialEnergyPivot(glp::Simplex &simplex);

    void chooseCol();
    void chooseRow();
};


#endif //GLPKTEST_POTENTIALENERGYPIVOT_H
