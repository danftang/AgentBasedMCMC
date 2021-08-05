//
// Created by daniel on 05/08/2021.
//

#ifndef GLPKTEST_CONVEXPMFPRODUCT_H
#define GLPKTEST_CONVEXPMFPRODUCT_H


#include <vector>
#include "ConvexPMF.h"

class ConvexPMFProduct {
    std::vector<ConvexPMF> atomicPMFs;

    operator ConvexPMF() {
        if(atomicPMFs.size() == 0) {
            return ConvexPMF([](const std::vector<double> &X) { return 0.0; }, 0);
        }
        return ConvexPMF([*this](const std::vector<double> &X) {
            double logP = 0.0;
            for(const ConvexPMF &pmf: atomicPMFs) logP += pmf.logP(X);
            return logP;
        },
                         atomicPMFs[0].nDimensions);
    }
};


#endif //GLPKTEST_CONVEXPMFPRODUCT_H
