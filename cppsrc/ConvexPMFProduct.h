//
// Created by daniel on 05/08/2021.
//

#ifndef GLPKTEST_CONVEXPMFPRODUCT_H
#define GLPKTEST_CONVEXPMFPRODUCT_H


#include <vector>
#include "ConvexPMF.h"

// If we're multiplying a very large number of terms, better to
// hold them separately to avoid stack overflow when optimisation
// is turned off or not possible.
//
// The cast operator allows this class to be cast to an ordinaary
// ConvexPMF when required.
class ConvexPMFProduct {
    std::vector<ConvexPMF> atomicPMFs;

    ConvexPMFProduct() { }

    ConvexPMFProduct(std::initializer_list<ConvexPMF> init): atomicPMFs(init) { }


    ConvexPMFProduct &operator *=(ConvexPMF atomicPMF) {
        atomicPMFs.push_back(std::move(atomicPMF));
        return *this;
    }


    ConvexPMFProduct &operator *=(ConvexPMFProduct others) {
        atomicPMFs.reserve(atomicPMFs.size() + others.atomicPMFs.size());
        for(ConvexPMF &atom: others.atomicPMFs) {
            atomicPMFs.push_back(std::move(atom));
        }
        return *this;
    }


    operator ConvexPMF() const & {
        if(atomicPMFs.size() == 0) return ConvexPMF::INVALID_PMF;
        return ConvexPMF([*this](const std::vector<double> &X) {
            double logP = 0.0;
            for(const ConvexPMF &pmf: atomicPMFs) logP += pmf.logP(X);
            return logP;
        },
                         atomicPMFs[0].nDimensions);
    }

    operator ConvexPMF() && {
        if(atomicPMFs.size() == 0) return ConvexPMF::INVALID_PMF;
        return ConvexPMF([atoms = std::move(atomicPMFs)](const std::vector<double> &X) {
            double logP = 0.0;
            for(const ConvexPMF &pmf: atoms) logP += pmf.logP(X);
            return logP;
            },
                         atomicPMFs[0].nDimensions);
    }

};



#endif //GLPKTEST_CONVEXPMFPRODUCT_H
