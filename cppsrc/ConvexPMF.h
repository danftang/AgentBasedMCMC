//
// Created by daniel on 26/07/2021.
//

#ifndef GLPKTEST_CONVEXPMF_H
#define GLPKTEST_CONVEXPMF_H

#include "glpkpp.h"
#include "PMF.h"
#include "SimplexMCMC.h"
#include "ConvexPolyhedron.h"

class ConvexPMF {
public:
    typedef SimplexMCMC DefaultSampler;

    ConvexPolyhedron            convexSupport; // all non-zero probability point lie on the vertices of this polyhedron
    std::function<double(const std::vector<double> &)>    logProb;   // function from vertex to log probability

    ConvexPMF(std::function<double(const std::vector<double> &)> logPrior): logProb(std::move(logPrior)) { }
    ConvexPMF(std::function<double(const std::vector<double> &)> logPrior,
              ConvexPolyhedron constraints): logProb(std::move(logPrior)), convexSupport(std::move(constraints)) { }


    double operator()(const std::vector<double> &X) const { return convexSupport.isValidSolution(X) ? logProb(X) : 0.0; }


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
