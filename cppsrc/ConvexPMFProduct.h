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
public:
    std::vector<ConvexPMF::PMF> atomicPMFs;
    ConvexPolyhedron support;
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

    ConvexPMFProduct &operator *=(ConvexPMF atom) {
        assert(nDimensions == atom.nDimensions);
        atomicPMFs.push_back(std::move(atom.logProb));
        support += std::move(atom.convexSupport);
        return *this;
    }


    ConvexPMFProduct &operator *=(ConvexPMFProduct others) {
        assert(nDimensions == others.nDimensions);
        atomicPMFs.reserve(atomicPMFs.size() + others.atomicPMFs.size());
        for(ConvexPMF::PMF &atom: others.atomicPMFs) {
            atomicPMFs.push_back(std::move(atom));
        }
        support += std::move(others.support);
        return *this;
    }

    ConvexPMFProduct operator *(ConvexPMF other) const & {
        assert(nDimensions == other.nDimensions);
        ConvexPMFProduct prod(nDimensions);
        prod.atomicPMFs.reserve(atomicPMFs.size()+1);
        for(const ConvexPMF::PMF &pmf: atomicPMFs) prod.atomicPMFs.push_back(pmf);
        prod.atomicPMFs.push_back(std::move(other.logProb));
        prod.support += std::move(other.convexSupport);
        return prod;
    }

    ConvexPMFProduct operator *(ConvexPMF other) && {
        (*this) *= other;
        return std::move(*this);
    }


    operator ConvexPMF() const & {
        return ConvexPMF([atoms = this->atomicPMFs](const std::vector<double> &X) {
            double logP = 0.0;
            for(const ConvexPMF::PMF &pmf: atoms) logP += pmf(X);
            return logP;
        }, nDimensions, support);
    }

    operator ConvexPMF() && {
        return ConvexPMF([atoms = std::move(atomicPMFs)](const std::vector<double> &X) {
            double logP = 0.0;
            for(const ConvexPMF::PMF &pmf: atoms) logP += pmf(X);
            return logP;
            }, nDimensions, std::move(support));
    }


};



#endif //GLPKTEST_CONVEXPMFPRODUCT_H
