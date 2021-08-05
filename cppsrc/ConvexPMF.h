//
// Created by daniel on 26/07/2021.
//

#ifndef GLPKTEST_CONVEXPMF_H
#define GLPKTEST_CONVEXPMF_H

#include "glpkpp.h"
#include "PMF.h"
#include "SimplexMCMC.h"
#include "ConvexPolyhedron.h"

// Represents a probability mass function defined on the vertices of a convex polyhedron.
class ConvexPMF {
public:
    typedef SimplexMCMC DefaultSampler; // put this in the Distribution class?

    std::function<double(const std::vector<double> &)>    logProb;   // function from vertex co-ord to log probability
    int                         nDimensions;    // the number of dimensions of the convex polyhedron of points
    ConvexPolyhedron            convexSupport;  // all non-zero probability points lie on the vertices of this polyhedron

//    ConvexPMF(std::function<double(const std::vector<double> &)> logPrior): logProb(std::move(logPrior)) { }
    ConvexPMF(
            std::function<double(const std::vector<double> &)> logPrior,
            int nDimensions,
            ConvexPolyhedron constraints = ConvexPolyhedron()
                      ):
              logProb(std::move(logPrior)),
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

//    virtual std::vector<double> nextSample() {
//        if(sampler == NULL) initSampler();
//        sampler->nextSample();
//        return sampler->X();
//    }


    ConvexPMF &operator *=(const ConvexPMF &other) {
        logProb = [logP = std::move(logProb), otherLogP = other.logProb](const std::vector<double> &X) {
            return logP(X) + otherLogP(X);
        };
        convexSupport += other.convexSupport;
        return *this;
    }


    ConvexPMF &operator *=(ConvexPMF &&other) {
        logProb = [logP = std::move(logProb), otherLogP = std::move(other.logProb)](const std::vector<double> &X) {
            return logP(X) + otherLogP(X);
        };
        convexSupport += std::move(other.convexSupport);
        return *this;
    }

    ConvexPMF operator *(const ConvexPMF &other) const & {
        ConvexPMF result(*this);
        result *= other;
        return result;
    }

    ConvexPMF operator *(const ConvexPMF &other) && {
        (*this) *= other;
        return std::move(*this);
    }


    // union of supports... removed because could be confusing to mix multiplication of PMF and addition of constraints(?)
//    ConvexPMF &operator *=(const std::vector<glp::Constraint> &otherSupport) {
//        convexSupport += otherSupport;
//        return *this;
//    }


protected:
//    SimplexMCMC *sampler = NULL;
//    glp::Problem *lpProblem = NULL;

//    void initSampler()          { sampler = new SimplexMCMC(convexSupport, logProb); }
//    void invalidateSampler()    {
//        if(sampler != NULL) {
//            delete sampler;
//            sampler = NULL;
//        }
//    }

};


#endif //GLPKTEST_CONVEXPMF_H
