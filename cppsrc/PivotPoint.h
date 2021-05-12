//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_PIVOTPOINT_H
#define GLPKTEST_PIVOTPOINT_H


#include <utility>
#include <tuple>

class PivotPoint {
public:
    int i;
    int j;
    std::vector<double> col;

    PivotPoint(int i, int j, std::vector<double> col): i(i), j(j), col(col) {}

//    operator std::tuple<int &,int &>() { return std::tie(i,j); }
};


#endif //GLPKTEST_PIVOTPOINT_H
