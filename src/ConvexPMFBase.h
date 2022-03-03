//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_CONVEXPMFBASE_H
#define GLPKTEST_CONVEXPMFBASE_H

#include <cmath>

template<typename T> class ConvexPMF;

template<typename DOMAIN>
class ConvexPMFBase {
protected:
    ConvexPMFBase(
        std::function<double(const DOMAIN &)> extendedLogP,
        int nDimensions,
        ConvexPolyhedron constraints
    ):
            extendedLogProb(std::move(extendedLogP)),
            nDimensions(nDimensions),
            convexSupport(std::move(constraints)) { }

public:
    typedef std::function<double(const DOMAIN &)> PMF;

    PMF                 extendedLogProb;// function from vertex co-ord to log probability, extended to give some "approximate" value in infeasible coords.
    int                 nDimensions;    // the number of dimensions of the convex polyhedron of points
    ConvexPolyhedron    convexSupport;  // all non-zero probability points lie on the vertices of this polyhedron

    double operator()(const DOMAIN &X) const { return logP(X); }

    bool isInSupport(const DOMAIN &X) const { return convexSupport.isValidSolution(X); }

    double logP(const DOMAIN &X) const {
        assert(X.size() == nDimensions);
        return isInSupport(X) ? extendedLogProb(X) : -INFINITY;
    }
    double P(const DOMAIN &X) const {
        assert(X.size() == nDimensions);
        return isInSupport(X) ? exp(extendedLogProb(X)) : 0.0;
    }

    // Multiplicatin of distributions. Equivalent to summation of logP
    // and union of constraints. If the PMFs have different nDimensions,
    // we implicitly increase the nDimensions of the lower dimensional distribution
    // to that of the higher by assuming the lower dimensional variables refer to the
    // prefix of the higher.
    ConvexPMF<DOMAIN> &operator *=(const ConvexPMF<DOMAIN> &other) {
        extendedLogProb = [extLogP = std::move(extendedLogProb), otherExtLogP = other.extendedLogProb](const DOMAIN &X) {
            return extLogP(X) + otherExtLogP(X);
        };
        convexSupport += other.convexSupport;
        nDimensions = std::max(nDimensions, other.nDimensions);
        return reinterpret_cast<ConvexPMF<DOMAIN> &>(*this);
    }


    ConvexPMF<DOMAIN> &operator *=(ConvexPMF<DOMAIN> &&other) {
        extendedLogProb = [extLogP = std::move(extendedLogProb), otherExtLogP = std::move(other.extendedLogProb)](const DOMAIN &X) {
            return extLogP(X) + otherExtLogP(X);
        };
        convexSupport += std::move(other.convexSupport);
        nDimensions = std::max(nDimensions, other.nDimensions);
        return reinterpret_cast<ConvexPMF<DOMAIN> &>(*this);
    }

    ConvexPMF<DOMAIN> operator *(const ConvexPMF<DOMAIN> &other) const & {
        ConvexPMF<DOMAIN> result(*this);
        result *= other;
        result.nDimensions = std::max(nDimensions, other.nDimensions);
        return result;
    }

    ConvexPMF<DOMAIN> operator *(const ConvexPMF<DOMAIN> &other) && {
        (*this) *= other;
        return std::move(*this);
    }
};


#endif //GLPKTEST_CONVEXPMFBASE_H
