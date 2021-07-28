//
// Created by daniel on 28/07/2021.
//

#ifndef GLPKTEST_DELTAPMF_H
#define GLPKTEST_DELTAPMF_H

#include "cmath"

class DeltaPMF {
public:
    std::vector<double> deltaPoint;

    DeltaPMF(std::vector<double> nonZeroPoint): deltaPoint(std::move(nonZeroPoint)) { }

    double logProb(const std::vector<double> &X) {
        return(std::equal(X.begin(), X.end(), deltaPoint.begin())?0:-INFINITY);
    }

    std::vector<double> nextSample() { return deltaPoint; }
};


#endif //GLPKTEST_DELTAPMF_H
