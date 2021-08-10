//
// Created by daniel on 06/08/2021.
//

#ifndef GLPKTEST_AGENTSTATEOBSERVATION_H
#define GLPKTEST_AGENTSTATEOBSERVATION_H

template<typename AGENT>
class AgentStateObservation {
public:
    State<AGENT> state;
    std::function<double(double)> logLikelihood; // probability of making the observation given the occupation number of 'state'
    double upperBound;
    double lowerBound;

//    AgentStateObservation(State<AGENT> state, ConvexPMF likelihood):
//    state(std::move(state)),
//    logLikelihood(std::move(likelihood)) {}

    AgentStateObservation(const State<AGENT> &state, int nObserved, double pObserveIfPresent):
    state(state),
    logLikelihood([p=pObserveIfPresent, m=nObserved](double n) { return log(boost::math::pdf(boost::math::binomial(n, p), m)); }),
    lowerBound(nObserved),
    upperBound(pObserveIfPresent==1.0?nObserved:INFINITY)
    {
    }

    double logP(double occupation) const {
        if(occupation < lowerBound || occupation > upperBound) return -INFINITY;
        return logLikelihood(occupation);
    }


    ConvexPolyhedron support() const {
        if(lowerBound > 0 || upperBound < state.occupationUpperBound()) {
            return ConvexPolyhedron({std::max(lowerBound,0.0) <= 1.0*state <= std::min(upperBound,state.occupationUpperBound())});
        }
        return ConvexPolyhedron();
    }
};


#endif //GLPKTEST_AGENTSTATEOBSERVATION_H
