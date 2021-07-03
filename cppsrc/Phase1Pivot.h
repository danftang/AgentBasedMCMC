//
// Created by daniel on 28/06/2021.
//

#ifndef GLPKTEST_PHASE1PIVOT_H
#define GLPKTEST_PHASE1PIVOT_H


#include "glpkpp.h"
#include "ProposalPivot.h"

class Phase1Pivot: public ProposalPivot {
public:
    explicit Phase1Pivot(glp::Simplex &simplex);
    Phase1Pivot(glp::Simplex &simplex, int i, int j): ProposalPivot(simplex, i, j) { }

    void chooseCol();
    void chooseRow();

    void setToPivotIndex(int pivotIndex);

protected:
    int setSimplexToInfeasibilityObjective();

    bool isActive(int pmfIndex);

    double infeasibility(double deltaj);

    double infeasibility();

    int infeasibilityCount();
};


#endif //GLPKTEST_PHASE1PIVOT_H
