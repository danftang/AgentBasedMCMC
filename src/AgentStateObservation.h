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
    ABM::occupation_type lowerBound;    // number actually observed
    double pObserveIfPresent;

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


    double logLikelihood(ABM::occupation_type realOccupation) const {
        if(realOccupation < lowerBound || realOccupation > upperBound()) return -std::numeric_limits<double>::infinity();
        return log(likelihood(realOccupation));
    }

    // likelihood function
    double likelihood(ABM::occupation_type realOccupation) const {
        return boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), lowerBound);
    }


    ABM::occupation_type upperBound() const { return pObserveIfPresent==1.0?lowerBound: state.fermionicOccupationUpperBound(); }


//    [[nodiscard]] ConvexPolyhedron<ABM::occupation_type> support() const {
//        if(lowerBound > 0 || upperBound() < state.occupationUpperBound()) {
//            return ConvexPolyhedron<ABM::occupation_type>({ lowerBound <= 1*state <= upperBound() });
//        }
//        return {};
//    }

    Constraint<ABM::occupation_type> constraint() const {
            return { lowerBound <= 1*state <= upperBound() };
    }

    std::function<double(ABM::occupation_type)> toLogProbFunction() const {
        return [copy = *this](ABM::occupation_type realOccupation) {
            return copy.logLikelihood(realOccupation);
        };
    }


    friend std::ostream &operator <<(std::ostream &out, const AgentStateObservation<AGENT> & observation) {
        out << observation.state << " n=" << observation.lowerBound << " p=" << observation.pObserveIfPresent;
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
