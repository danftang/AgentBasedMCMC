//
// Created by daniel on 17/05/2021.
//
// Holds calculated information on a column during pivot operations

#ifndef GLPKTEST_PHASE2PIVOT_H
#define GLPKTEST_PHASE2PIVOT_H


#include <vector>
#include "ProposalPivot.h"
#include "SimplexMCMC.h"

class Phase2Pivot: public ProposalPivot {
public:
//    static constexpr double tol = SimplexMCMC::tol; // tolerance to consider a value zero

    std::vector<int>    pivotRows;  // rows on which this column can be pivoted while maintaining infeasibility (structural vars precede auxiliary)
    int                 nStructuralPivotRows; // Number of entries in pivotRows that correspond to structural pivots (i.e. not auxiliary)

    Phase2Pivot(glp::Simplex &simplex, int j): Phase2Pivot(simplex, j, simplex.tableauCol(j)) { }
    Phase2Pivot(glp::Simplex &simplex, int j, std::vector<double> column); // j = column to pivot on


    Phase2Pivot reverse(glp::Simplex &lp) const;
    bool isDegenerate() const { return fabs(deltaj) < tol; }

protected:
    void orderPivotRows(glp::Simplex &lp);
};


#endif //GLPKTEST_PHASE2PIVOT_H
