//
// Created by daniel on 26/07/2021.
//

#ifndef GLPKTEST_CONVEXPMF_H
#define GLPKTEST_CONVEXPMF_H

#include "glpkpp.h"
#include "PMF.h"
#include "SimplexMCMC.h"

class ConvexPMF {
public:
    glp::Problem    convexSupport; // all non-zero probability point lie on the vertices of this polyhedron
    std::function<double(const std::vector<double> &)>    logPrior;   // function from vertex to log probability
    SimplexMCMC *sampler;


    ~ConvexPMF() {
        if(sampler != NULL) delete sampler;
    }


    double logProb(const std::vector<double> &X) const { return convexSupport.isValidSolution(X)?logPrior(X):0.0; }


    std::vector<double> nextSample() {
        if(sampler == NULL) initSampler();
        sampler->nextSample();
        return sampler->X();
    }


    ConvexPMF &operator *=(const ConvexPMF &other) {
        logPrior = [logP = std::move(logPrior), otherLogP = other.logPrior](const std::vector<double> &X) {
            return logP(X) + otherLogP(X);
        };
        (*this) *= other.convexSupport;
        return *this;
    }


    ConvexPMF &operator *=(ConvexPMF &&other) {
        logPrior = [logP = std::move(logPrior), otherLogP = std::move(other.logPrior)](const std::vector<double> &X) {
            return logP(X) + otherLogP(X);
        };
        (*this) *= other.convexSupport;
        return *this;
    }


    // union of supports
    ConvexPMF &operator *=(const glp::Problem &otherSupport) {
        for(int c=1; c<=otherSupport.nConstraints(); ++c) {
            convexSupport.addConstraint(otherSupport.getConstraint(c));
        }
        return *this;
    }


protected:
    void initSampler()          { sampler = new SimplexMCMC(convexSupport, logPrior); }
    void invalidateSampler()    { delete sampler; sampler = NULL; }

};


#endif //GLPKTEST_CONVEXPMF_H
