//
// Created by daniel on 17/05/2021.
//

#ifndef GLPKTEST_COLUMNPIVOT_H
#define GLPKTEST_COLUMNPIVOT_H


#include <vector>
#include "Pivot.h"

class ColumnPivot: public Pivot {
public:
    std::vector<int>    pivotRows;

    ColumnPivot(glp::Simplex &simplex, int j): ColumnPivot(simplex, j, simplex.tableauCol(j)) { }
    ColumnPivot(glp::Simplex &simplex, int j, std::vector<double> column); // j = column to pivot on

    ColumnPivot reverse(glp::Simplex &lp) const;
};


#endif //GLPKTEST_COLUMNPIVOT_H
