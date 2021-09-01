//
// Created by daniel on 06/08/2021.
//

#ifndef GLPKTEST_AGENTSTATEOBSERVATION_H
#define GLPKTEST_AGENTSTATEOBSERVATION_H

#include "boost/math/distributions/binomial.hpp"
#include "ConvexPMFProduct.h"
#include "TrajectoryPriorDistribution.h"
#include "debug.h"

template<typename AGENT>
class AgentStateObservation {
public:
    State<AGENT> state;
//    std::function<double(double)> logLikelihood; // probability of making the observation given the occupation number of 'state'
//    double upperBound;
    double lowerBound;
    double pObserveIfPresent;
    double logPExpectation;
    double normalisation;

//    AgentStateObservation(State<AGENT> state, ConvexPMF likelihood):
//    state(std::move(state)),
//    logLikelihood(std::move(likelihood)) {}

    AgentStateObservation(const State<AGENT> &state, int nObserved, double pObserveIfPresent):
    state(state),
//    logLikelihood([p=pObserveIfPresent, m=nObserved](double n) { return log(boost::math::pdf(boost::math::binomial(n, p), m)); }),
    lowerBound(nObserved),
//    upperBound(pObserveIfPresent==1.0?nObserved:INFINITY),
    pObserveIfPresent(pObserveIfPresent),
    normalisation(1.0)
    {
        double p;
//        double sumPsq = 0.0;
        double sumP = 0.0;
        int n=nObserved;
        do {
//            std::cout << "P(" << n << ") = ";
            p = P(n++);
//            std::cout << p << std::endl;
//            sumPsq += p*p;
            sumP += p;
        } while(p > 1e-6);
        normalisation = 1.0/sumP;
        logPExpectation = log(P(nObserved+1.0));
        debug(std::cout << "Generating observation " << state << " " << nObserved << " " << pObserveIfPresent << " " << exp(logPExpectation) << std::endl);
    }

    double extendedLogP(double realOccupation) const {
        double rRealOccupation = std::round(realOccupation);
        if(rRealOccupation < lowerBound || rRealOccupation > upperBound()) {
            return logPExpectation; // case when realOccupation == lowerBound TODO: make a better approximation to this?
        }
        return log(P(rRealOccupation));
    }

    double P(double realOccupation) const {
        return normalisation*boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), lowerBound);
    }


    double upperBound() const { return pObserveIfPresent==1.0?lowerBound:INFINITY; }


    ConvexPolyhedron support() const {
        if(lowerBound > 0 || upperBound() < state.occupationUpperBound()) {
            return ConvexPolyhedron({std::max(lowerBound,0.0) <= 1.0*state <= std::min(upperBound(),state.occupationUpperBound())});
        }
        return ConvexPolyhedron();
    }



};


#endif //GLPKTEST_AGENTSTATEOBSERVATION_H
