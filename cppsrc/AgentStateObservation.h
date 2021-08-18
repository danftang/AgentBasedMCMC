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

//    AgentStateObservation(State<AGENT> state, ConvexPMF likelihood):
//    state(std::move(state)),
//    logLikelihood(std::move(likelihood)) {}

    AgentStateObservation(const State<AGENT> &state, int nObserved, double pObserveIfPresent):
    state(state),
//    logLikelihood([p=pObserveIfPresent, m=nObserved](double n) { return log(boost::math::pdf(boost::math::binomial(n, p), m)); }),
    lowerBound(nObserved),
//    upperBound(pObserveIfPresent==1.0?nObserved:INFINITY),
    pObserveIfPresent(pObserveIfPresent)
    {
        debug(std::cout << "Generating observation " << state << " " << nObserved << " " << pObserveIfPresent << std::endl);
    }

    double logP(double realOccupation) const {
        if(realOccupation < lowerBound || realOccupation > upperBound()) return -INFINITY;
        return log(boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), lowerBound));// logLikelihood(occupation);
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
