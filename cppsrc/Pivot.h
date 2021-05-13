//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_PIVOT_H
#define GLPKTEST_PIVOT_H


#include <utility>
#include <tuple>
#include <vector>

#include "glpkpp.h"

class Pivot {
public:

    int i;
    int j;
    std::vector<double> col;
    std::vector<int>    pivotRows;

//    Pivot(int i, int j, std::vector<double> col, std::vector<int> pivotRows): i(i), j(j), col(std::move(col)), pivotRows(std::move(pivotRows)) {}

    Pivot() = default;
    Pivot(glp::Simplex &simplex, int j); // j = column to pivot on
    Pivot(glp::Simplex &lp, int i, int j, std::vector<double> col);

    //    operator std::tuple<int &,int &>() { return std::tie(i,j); }
    std::vector<double> reverseCol();
};



#endif //GLPKTEST_PIVOT_H
