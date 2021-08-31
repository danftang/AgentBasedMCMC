//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_CONVEXPMFBASE_H
#define GLPKTEST_CONVEXPMFBASE_H

template<typename T> class ConvexPMF;

template<typename DOMAIN>
class ConvexPMFBase {
protected:
    ConvexPMFBase(
        std::function<double(const DOMAIN &)> logP,
        int nDimensions,
        ConvexPolyhedron constraints
    ):
    logProb(std::move(logP)),
    nDimensions(nDimensions),
    convexSupport(std::move(constraints)) { }

public:
    typedef std::function<double(const DOMAIN &)> PMF;

    PMF                 logProb;        // function from vertex co-ord to log probability
    int                 nDimensions;    // the number of dimensions of the convex polyhedron of points
    ConvexPolyhedron    convexSupport;  // all non-zero probability points lie on the vertices of this polyhedron

    double operator()(const DOMAIN &X) const { return logP(X); }

    bool isInSupport(const DOMAIN &X) const { return convexSupport.isValidSolution(X); }

    double logP(const DOMAIN &X) const {
        assert(X.size() == nDimensions);
        return isInSupport(X) ? logProb(X) : 0.0;
    }
    double P(const DOMAIN &X) const { return exp(logP(X)); }

    // Multiplicatin of distributions. Equivalent to summation of logP
    // and union of constraints. If the PMFs have different nDimensions,
    // we implicitly increase the nDimensions of the lower dimensional distribution
    // to that of the higher by assuming the lower dimensional variables refer to the
    // prefix of the higher.
    ConvexPMF<DOMAIN> &operator *=(const ConvexPMF<DOMAIN> &other) {
        logProb = [logP = std::move(logProb), otherLogP = other.logProb](const DOMAIN &X) {
            return logP(X) + otherLogP(X);
        };
        convexSupport += other.convexSupport;
        nDimensions = std::max(nDimensions, other.nDimensions);
        return reinterpret_cast<ConvexPMF<DOMAIN> &>(*this);
    }


    ConvexPMF<DOMAIN> &operator *=(ConvexPMF<DOMAIN> &&other) {
        logProb = [logP = std::move(logProb), otherLogP = std::move(other.logProb)](const DOMAIN &X) {
            return logP(X) + otherLogP(X);
        };
        convexSupport += std::move(other.convexSupport);
        nDimensions = std::max(nDimensions, other.nDimensions);
        return *this;
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
