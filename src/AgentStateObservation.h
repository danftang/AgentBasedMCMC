//
// Created by daniel on 06/08/2021.
//

#ifndef GLPKTEST_AGENTSTATEOBSERVATION_H
#define GLPKTEST_AGENTSTATEOBSERVATION_H

#include "boost/math/distributions/binomial.hpp"
#include "ConvexPolyhedron.h"
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

    AgentStateObservation() {};

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
        debug(std::cout << "Generating observation " << *this << std::endl);
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

    static std::vector<AgentStateObservation<AGENT>>
    generateObservations(const Trajectory<AGENT> &realTrajectory, double pMakeObservation, double pObserveIfPresent) {
        int nTimesteps = realTrajectory.nTimesteps();
        std::vector<AgentStateObservation<AGENT>> observations;
        for (int t=0; t<nTimesteps;++t) {
            for (int agentId=0; agentId < AGENT::domainSize(); ++agentId) {
                if (Random::nextDouble() < pMakeObservation) {
                    AGENT agent(agentId);
                    int nObserved = Random::nextBinomial(realTrajectory(t,agent), pObserveIfPresent);
                    observations.emplace_back(State<AGENT>(t,agent), nObserved, pObserveIfPresent);
                    assert(observations.back().support().isValidSolution(realTrajectory));
                }
            }
        }
        return observations;
    }

    friend std::ostream &operator <<(std::ostream &out, const AgentStateObservation<AGENT> & observation) {
        out << observation.state << " " << observation.lowerBound << " " << observation.pObserveIfPresent << " " << exp(observation.logPExpectation);
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & state & lowerBound & pObserveIfPresent & logPExpectation & normalisation;
    }


};


#endif //GLPKTEST_AGENTSTATEOBSERVATION_H
