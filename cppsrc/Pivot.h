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
//    std::vector<int>    pivotRows;

    Pivot() = default;
//    Pivot(glp::Simplex &simplex, int j); // j = column to pivot on
    Pivot(int i, int j, std::vector<double> column): i(i), j(j), col(std::move(column)) { }


    //    operator std::tuple<int &,int &>() { return std::tie(i,j); }
    std::vector<double> reverseCol() const;
//    Pivot reverse() const { return Pivot(i, j, reverseCol()); }
};



#endif //GLPKTEST_PIVOT_H
