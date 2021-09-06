//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_BERNOULLIMODELSTATE_H
#define GLPKTEST_BERNOULLIMODELSTATE_H

#include "Distribution.h"

template<typename AGENT>
class BernoulliModelState: public Distribution<ModelState<AGENT>> {
public:
//    static constexpr double infeasibleP = 0.5;
    std::function<double(AGENT)> prob;

    BernoulliModelState(std::function<double(AGENT)> probabilities): prob(probabilities) {
    }

    BernoulliModelState(std::vector<double> probabilities): prob([p = std::move(probabilities)](AGENT agent) {
        assert(agent < p.size());
        return p[agent];
    }) { }


    ConvexPMF<ModelState<AGENT>> PMF() const {
        return ConvexPMF<ModelState<AGENT>>(
                [*this](const ModelState<AGENT> &M) {
                    double extendedLogP = 0.0;
                    for(int agentId = 0; agentId < AGENT::domainSize(); ++agentId) {
                        int occupation = std::round(M[agentId]);
                        double p = prob(AGENT(agentId));
                        switch(occupation) {
                            case 0: extendedLogP += log(1.0 - p); break;
                            case 1: extendedLogP += log(p); break;
                            default: extendedLogP += log(infeasibleExpectationFraction *(p*p + (1.0-p)*(1.0-p)));
                        }
                    }
                    return extendedLogP;
                },
                contraints()
                );
    }

    std::function<ModelState<AGENT>()> sampler() const {
        return [*this]() { return nextSample(); };
    }

    ModelState<AGENT> nextSample() const {
        ModelState<AGENT> sample;
        for(int agentId = 0; agentId < nDimensions(); ++agentId) {
            sample[agentId] = Random::nextDouble() < prob(agentId)?1.0:0.0;
        }
        return sample;
    }

    ConvexPolyhedron contraints() const {
        ConvexPolyhedron constraints;
        constraints.reserve(nDimensions());
        for(int agentId=0; agentId<nDimensions(); ++agentId) {
            double p = prob(agentId);
            double lowerBound = p==1.0?1.0:0.0;
            double upperBound = p==0.0?0.0:1.0;
            constraints.push_back(lowerBound <= 1.0*glp::X(agentId) <= upperBound);
        }
        return constraints;
    }

    int nDimensions() const { return AGENT::domainSize(); }

    friend std::ostream &operator <<(std::ostream &out, const BernoulliModelState<AGENT> &pmf) {
        out << "{ ";
        for(int agentId=0; agentId < pmf.nDimensions(); ++agentId) {
            out << pmf.prob(agentId) << " ";
        }
        out << "}";
        return out;
    }
};


#endif //GLPKTEST_BERNOULLIMODELSTATE_H
