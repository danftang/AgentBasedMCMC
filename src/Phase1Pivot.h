//
// Created by daniel on 28/06/2021.
//

#ifndef GLPKTEST_PHASE1PIVOT_H
#define GLPKTEST_PHASE1PIVOT_H

#include "ProposalPivot.h"

class Phase1Pivot: public ProposalPivot {
public:
    explicit Phase1Pivot(SimplexMCMC &simplex);
    Phase1Pivot(SimplexMCMC &simplex, int i, int j);

    std::vector<double> infeasibilityGradient;  // infeasibility objective by row [1...nBasic()]

    void chooseCol();
    void chooseRow();

//    void setToPivotIndex(int pivotIndex);

protected:
//    int initInfeasibilityGradient();

//    bool isActive(int pmfIndex);

    double infeasibility(double deltaj);

    double infeasibility();

    int infeasibilityCount();
};


#endif //GLPKTEST_PHASE1PIVOT_H
