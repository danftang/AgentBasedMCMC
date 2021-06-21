//
// Created by daniel on 17/05/2021.
//
// Holds calculated information on a column during pivot operations

#ifndef GLPKTEST_COLUMNPIVOT_H
#define GLPKTEST_COLUMNPIVOT_H


#include <vector>
#include "ProposalPivot.h"

class ColumnPivot: public ProposalPivot {
public:
    static constexpr double tol = 1e-8; // tolerance to consider a value zero

    std::vector<int>    pivotRows;  // rows on which this column can be pivoted while maintaining infeasibility (structural vars precede auxiliary)
    int                 nStructuralPivotRows; // Number of entries in pivotRows that correspond to structural pivots (i.e. not auxiliary)

    ColumnPivot(glp::Simplex &simplex, int j): ColumnPivot(simplex, j, simplex.tableauCol(j)) { }
    ColumnPivot(glp::Simplex &simplex, int j, std::vector<double> column); // j = column to pivot on


    ColumnPivot reverse(glp::Simplex &lp) const;
    bool isDegenerate() const { return fabs(deltaj) < tol; }

protected:
    void orderPivotRows(glp::Simplex &lp);
};


#endif //GLPKTEST_COLUMNPIVOT_H
