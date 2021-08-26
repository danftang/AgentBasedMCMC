//
// Created by daniel on 26/07/2021.
//

#ifndef GLPKTEST_CONVEXPMF_H
#define GLPKTEST_CONVEXPMF_H

#include "glpkpp.h"
#include "SimplexMCMC.h"
#include "ConvexPolyhedron.h"

// Represents a probability mass function defined on the vertices of a convex polyhedron.
class ConvexPMF {
public:
    typedef SimplexMCMC DefaultSampler; // put this in the Distribution class?
    typedef std::function<double(const std::vector<double> &)> PMF;

    static const ConvexPMF INVALID_PMF;

    PMF                 logProb;        // function from vertex co-ord to log probability
    int                 nDimensions;    // the number of dimensions of the convex polyhedron of points
    ConvexPolyhedron    convexSupport;  // all non-zero probability points lie on the vertices of this polyhedron

    ConvexPMF(
            PMF logP,
            int nDimensions,
            ConvexPolyhedron constraints = ConvexPolyhedron()
                      ):
              logProb(std::move(logP)),
              nDimensions(nDimensions),
              convexSupport(std::move(constraints)) { }


    double operator()(const std::vector<double> &X) const { return convexSupport.isValidSolution(X) ? logProb(X) : 0.0; }

    double logP(const std::vector<double> &X) const {
        assert(X.size() == nDimensions);
        return convexSupport.isValidSolution(X) ? logProb(X) : 0.0;
    }
    double P(const std::vector<double> &X) const { return exp(logP(X)); }

    std::function<std::vector<double>()> sampler() {
        return [mcmc = SimplexMCMC(*this)]() mutable { return std::vector<double>(mcmc.nextSample()); };
    }


    // Multiplicatin of distributions. Equivalent to summation of logP
    // and union of constraints. If the PMFs have different nDimensions,
    // we implicitly increase the nDimensions of the lower dimensional distribution
    // to that of the higher by assuming the lower dimensional variables refer to the
    // prefix of the higher.
    ConvexPMF &operator *=(const ConvexPMF &other) {
        logProb = [logP = std::move(logProb), otherLogP = other.logProb](const std::vector<double> &X) {
            return logP(X) + otherLogP(X);
        };
        convexSupport += other.convexSupport;
        nDimensions = std::max(nDimensions, other.nDimensions);
        return *this;
    }


    ConvexPMF &operator *=(ConvexPMF &&other) {
        logProb = [logP = std::move(logProb), otherLogP = std::move(other.logProb)](const std::vector<double> &X) {
            return logP(X) + otherLogP(X);
        };
        convexSupport += std::move(other.convexSupport);
        nDimensions = std::max(nDimensions, other.nDimensions);
        return *this;
    }

    ConvexPMF operator *(const ConvexPMF &other) const & {
        ConvexPMF result(*this);
        result *= other;
        result.nDimensions = std::max(nDimensions, other.nDimensions);
        return result;
    }

    ConvexPMF operator *(const ConvexPMF &other) && {
        (*this) *= other;
        return std::move(*this);
    }

//    // Standard 1-dimensional PMFs
//
//    static ConvexPMF binomial(int nTrials, double p);
//    static ConvexPMF delta(int n);



protected:

};


#endif //GLPKTEST_CONVEXPMF_H
