//
// Created by daniel on 28/07/2021.
//

#ifndef GLPKTEST_DELTAPMF_H
#define GLPKTEST_DELTAPMF_H

#include "cmath"

class DeltaPMF: public ConvexPMF {
public:
    typedef DeltaPMF DefaultSampler; // copy or move construction as a sampler

    std::vector<double> deltaPoint;

    DeltaPMF(std::vector<double> nonZeroPoint):
    ConvexPMF([this](const std::vector<double> &X) { return (*this)(X); }),
    deltaPoint(std::move(nonZeroPoint)) {
        convexSupport.reserve(deltaPoint.size()-1);
        for(int i=0; i<deltaPoint.size(); ++i) {
            convexSupport.push_back(1.0*glp::X(i) == deltaPoint[i]);
        }
    }

    double operator()(const std::vector<double> &X) {
        return(std::equal(X.begin(), X.end(), deltaPoint.begin())?0.0:-INFINITY);
    }

    std::vector<double> nextSample() { return deltaPoint; }

};


#endif //GLPKTEST_DELTAPMF_H
