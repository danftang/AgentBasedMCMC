//
// Created by daniel on 11/11/2021.
//

#ifndef GLPKTEST_POISSONMODELSTATE_H
#define GLPKTEST_POISSONMODELSTATE_H

#include "Distribution.h"

template<typename AGENT>
class PoissonModelState: public Distribution<ModelState<AGENT>> {
public:
    std::function<double(AGENT)> poissonLambda; // function from agent state to Poisson lambda for that agent

    PoissonModelState(std::function<double(AGENT)> poissonLambda): poissonLambda(poissonLambda) {
    }

    PoissonModelState(std::vector<double> poissonLambdas): poissonLambda([p = std::move(poissonLambdas)](AGENT agent) {
        assert(agent < p.size());
        return p[agent];
    }) { }


    ConvexPMF<ModelState<AGENT>> PMF() const {
        return ConvexPMF<ModelState<AGENT>>(
                [*this](const ModelState<AGENT> &M) {
                    double logP = 0.0;
                    for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                        double k = std::max(0.0,std::round(M[agentId]));
                        double l = poissonLambda(AGENT(agentId));
                        logP += k*log(l) - l - lgamma(k+1); // log of Poisson
                    }
                    return logP;
                },
                contraints()
        );
    }

    std::function<ModelState<AGENT>()> sampler() const {
        return [*this]() { return nextSample(); };
    }

    ModelState<AGENT> nextSample() const {
        ModelState<AGENT> state;
        for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
            state[agentId] = Random::nextPoisson(poissonLambda(AGENT(agentId)));
        }
        return state;
    }

    ConvexPolyhedron contraints() const {
        return ConvexPolyhedron();
    }

    int nDimensions() const { return AGENT::domainSize(); }

    friend std::ostream &operator<<(std::ostream &out, const PoissonModelState <AGENT> &pmf) {
        out << "{ ";
        for (int agentId = 0; agentId < pmf.nDimensions(); ++agentId) {
            out << pmf.poissonLambda(AGENT(agentId)) << " ";
        }
        out << "}";
        return out;
    }

};


#endif //GLPKTEST_POISSONMODELSTATE_H
