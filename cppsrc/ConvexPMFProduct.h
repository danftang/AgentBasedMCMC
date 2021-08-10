//
// Created by daniel on 05/08/2021.
//

#ifndef GLPKTEST_CONVEXPMFPRODUCT_H
#define GLPKTEST_CONVEXPMFPRODUCT_H


#include <vector>
#include "ConvexPMF.h"
#include "PMFProduct.h"

// If we're multiplying a very large number of terms, better to
// hold them separately to avoid stack overflow when optimisation
// is turned off or not possible.
//
// The cast operator allows this class to be cast to an ordinaary
// ConvexPMF when required.
class ConvexPMFProduct {
protected:
    PMFProduct atomicPMFs;
public:
    ConvexPolyhedron convexSupport;
    int nDimensions;

    ConvexPMFProduct(int nDimensions): nDimensions(nDimensions) { }

//    ConvexPMFProduct(std::initializer_list<ConvexPMF> init): atomicPMFs(init.size()) {
//        assert(init.size() != 0);
//        nDimensions = init.begin()->nDimensions;
//    }
//
//    ConvexPMFProduct(std::vector<ConvexPMF> atoms): atomicPMFs(std::move(atoms)) {
//        assert(atoms.size() != 0);
//        nDimensions = atomicPMFs[0].nDimensions;
//    }

    int size() const { return atomicPMFs.size(); }

    double logP(const std::vector<double> &X) const {
        return atomicPMFs.logP(X);
    }

    ConvexPMFProduct &operator *=(ConvexPMF atom) {
        assert(nDimensions == atom.nDimensions);
        atomicPMFs *= std::move(atom.logProb);
        convexSupport += std::move(atom.convexSupport);
        return *this;
    }


    ConvexPMFProduct &operator *=(ConvexPMFProduct others) {
        assert(nDimensions == others.nDimensions);
        atomicPMFs *= others.atomicPMFs;
        convexSupport += std::move(others.convexSupport);
        return *this;
    }

    ConvexPMFProduct operator *(ConvexPMF other) const & {
        assert(nDimensions == other.nDimensions);
        ConvexPMFProduct prod(nDimensions);
        prod.atomicPMFs = atomicPMFs * std::move(other.logProb);
        prod.convexSupport.reserve(convexSupport.size() + other.convexSupport.size());
        prod.convexSupport += convexSupport;
        prod.convexSupport += std::move(other.convexSupport);
        return prod;
    }

    ConvexPMFProduct operator *(ConvexPMF other) && {
        (*this) *= other;
        return std::move(*this);
    }


    operator ConvexPMF() const & { return ConvexPMF(atomicPMFs.PMF(), nDimensions, convexSupport); }

    operator ConvexPMF() && {
        return ConvexPMF(std::move(atomicPMFs).PMF(), nDimensions, std::move(convexSupport));
    }


};



#endif //GLPKTEST_CONVEXPMFPRODUCT_H
