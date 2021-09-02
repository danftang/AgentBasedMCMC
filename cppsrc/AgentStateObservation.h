//
// Created by daniel on 06/08/2021.
//

#ifndef GLPKTEST_AGENTSTATEOBSERVATION_H
#define GLPKTEST_AGENTSTATEOBSERVATION_H

#include "boost/math/distributions/binomial.hpp"
#include "ConvexPMFProduct.h"
#include "TrajectoryPriorDistribution.h"
#include "debug.h"
#include "constants.h"

template<typename AGENT>
class AgentStateObservation {
public:
    State<AGENT> state;
    double lowerBound;
    double pObserveIfPresent;
    double logPExpectation;
    double normalisation;


    AgentStateObservation(const State<AGENT> &state, int nObserved, double pObserveIfPresent):
    state(state),
    lowerBound(nObserved),
    pObserveIfPresent(pObserveIfPresent),
    normalisation(1.0)
    {
        double p;
        double sumPsq = 0.0;
        double sumP = 0.0;
        int n=nObserved;
        do {
//            std::cout << "P(" << n << ") = ";
            p = P(n++);
//            std::cout << p << std::endl;
            sumPsq += p*p;
            sumP += p;
        } while(p > 1e-6);
        logPExpectation = log(infeasibleExpectationFraction*sumPsq/(sumP*sumP));
        normalisation = 1.0/sumP;
//        logPExpectation = log(P(nObserved+1.0));
        debug(std::cout << "Generating observation " << state << " " << nObserved << " " << pObserveIfPresent << " " << exp(logPExpectation) << std::endl);
    }


    double extendedLogP(double realOccupation) const {
        double rRealOccupation = std::round(realOccupation);
        if(rRealOccupation < lowerBound || rRealOccupation > upperBound()) return logPExpectation;
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
