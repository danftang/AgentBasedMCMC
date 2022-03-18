// Likelihood function for an observation of the number of agents in a given state
// at a given time.
// The observation can be made imperfect in that present agents may not be
// observed, but observed agents are always present (i.e. no false positives).
//
// So, the likelihood function is the Binomial distribution, with a fixed number of
// observed events, but variable number of trials
//
// Created by daniel on 06/08/2021.
//

#ifndef GLPKTEST_AGENTSTATEOBSERVATION_H
#define GLPKTEST_AGENTSTATEOBSERVATION_H

#include "boost/math/distributions/binomial.hpp"
#include "ConvexPolyhedron.h"
#include "debug.h"
#include "ABM.h"

template<typename AGENT>
class AgentStateObservation {
public:
    State<AGENT> state;         // state that was observed
    const double lowerBound;    // number actually observed
    const double pObserveIfPresent;

    AgentStateObservation() {};

    AgentStateObservation(const State<AGENT> &state, int nObserved, double pObserveIfPresent):
    state(state),
    lowerBound(nObserved),
    pObserveIfPresent(pObserveIfPresent)
    {
        assert(nObserved >= 0);
        assert(pObserveIfPresent > 0.0 && pObserveIfPresent <= 1.0);
//        double p;
//        double sumPsq = 0.0;
//        double sumP = 0.0;
//        int n=nObserved;
//        do {
//            p = P(n++);
//            sumPsq += p*p;
//            sumP += p;
//        } while(p > 1e-6);
//        logPExpectation = log(infeasibleExpectationFraction*sumPsq/(sumP*sumP));
        debug(std::cout << "Generating observation " << *this << std::endl);
    }


    double logLikelihood(double realOccupation) const {
        double rRealOccupation = std::round(realOccupation);
        if(rRealOccupation < lowerBound || rRealOccupation > upperBound()) return -INFINITY;
        return log(likelihood(rRealOccupation));
    }

    // likelihood function
    double likelihood(double realOccupation) const {
        return boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), lowerBound);
    }


    double upperBound() const { return pObserveIfPresent==1.0?lowerBound:state.occupationUpperBound(); }


    [[nodiscard]] ConvexPolyhedron<ABM::occupation_type> support() const {
        if(lowerBound > 0 || upperBound() < state.occupationUpperBound()) {
            return ConvexPolyhedron<ABM::occupation_type>({ lowerBound <= 1*state <= upperBound() });
        }
        return {};
    }


    friend std::ostream &operator <<(std::ostream &out, const AgentStateObservation<AGENT> & observation) {
        out << observation.state << " " << observation.lowerBound << " " << observation.pObserveIfPresent << " " << exp(observation.logPExpectation);
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & state & lowerBound & pObserveIfPresent;
    }


};


#endif //GLPKTEST_AGENTSTATEOBSERVATION_H
