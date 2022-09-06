// Likelihood function for an observation of the number of agents in a given state
// at a given time.
// The observation can be made noisy in that present agents may not be
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
#include "ConstrainedFactorisedDistribution.h"

template<typename AGENT>
class NoisyAgentStateObservation {
public:
    const State<AGENT> state;                 // state that was observed
    const ABM::occupation_type nObserved;     // number actually observed
    const double pObserveIfPresent;           // probability of observing an agent if it is present

    NoisyAgentStateObservation() {};

    NoisyAgentStateObservation(const State<AGENT> &state, ABM::occupation_type nObserved, double pObserveIfPresent):
            state(state),
            nObserved(nObserved),
            pObserveIfPresent(pObserveIfPresent)
    {
        assert(nObserved >= 0);
        assert(pObserveIfPresent > 0.0 && pObserveIfPresent < 1.0); // if pObserveIfPresent = 1.0 use NoiselessAgentStateObservation
//        debug(std::cout << "Generating observation " << *this << std::endl);
    }

    SparseWidenedFunction<double, const Trajectory<AGENT> &> toSparseWidenedFunction() {
        return SparseWidenedFunction<double, const Trajectory<AGENT> &>(
        [state = state,nObserved = nObserved,pObserveIfPresent = pObserveIfPresent](const Trajectory<AGENT> &trajectory) {
            return widenedLogLikelihood(state, nObserved, pObserveIfPresent, trajectory);
        },
        state.forwardOccupationDependencies()
        );
    }

    static std::pair<double,bool> widenedLogLikelihood(const State<AGENT> &state, ABM::occupation_type nObserved, double pObserveIfPresent, const Trajectory<AGENT> &trajecotry) {
        ABM::occupation_type realOccupation = trajecotry[state];
        if(realOccupation < nObserved) {
            return std::pair(nObserved*log(pObserveIfPresent) - ABM::kappa*(nObserved - realOccupation),false);
        }
        return std::pair(log(boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), nObserved)),true);
    }



//    double logLikelihood(ABM::occupation_type realOccupation) const {
//        if(realOccupation < lowerBound || realOccupation > upperBound()) return -std::numeric_limits<double>::infinity();
//        return log(likelihood(realOccupation));
//    }
//
//    // likelihood function
//    double likelihood(ABM::occupation_type realOccupation) const {
//        return boost::math::pdf(boost::math::binomial(realOccupation, pObserveIfPresent), lowerBound);
//    }
//
//
//    ABM::occupation_type upperBound() const { return pObserveIfPresent==1.0?lowerBound: state.fermionicOccupationUpperBound(); }
//
//
//    Constraint<ABM::occupation_type> constraint() const {
//            return { lowerBound <= 1*state <= upperBound() };
//    }
//
//    std::function<double(ABM::occupation_type)> toLogProbFunction() const {
//        return [copy = *this](ABM::occupation_type realOccupation) {
//            return copy.logLikelihood(realOccupation);
//        };
//    }


    friend std::ostream &operator <<(std::ostream &out, const NoisyAgentStateObservation<AGENT> & observation) {
        out << observation.state << " n=" << observation.nObserved << " p=" << observation.pObserveIfPresent;
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & state & nObserved & pObserveIfPresent;
    }


};


#endif //GLPKTEST_AGENTSTATEOBSERVATION_H
