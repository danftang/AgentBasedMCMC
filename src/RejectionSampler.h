//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H

#include <functional>
#include "ModelState.h"
#include "Prior.h"
#include "Likelihood.h"

template<class DOMAIN>
class RejectionSampler {
public:
    std::function<DOMAIN()>                 priorSampler;
    std::function<double(const DOMAIN &)>   likelihood;

    RejectionSampler(std::function<DOMAIN()> PriorSampler, std::function<double(const DOMAIN &)> Likelihood):
        priorSampler(std::move(PriorSampler)),
        likelihood(std::move(Likelihood)) {
    }


    template<class STARTSTATE, class AGENT>
    RejectionSampler(Prior<STARTSTATE> & Prior, Likelihood<AGENT> & Likelihood) :
        priorSampler([&Prior]() { return Prior.nextSample(); }),
        likelihood([&Likelihood](const Trajectory<AGENT> &X) {
            double logP = Likelihood.logPexact(X);
            return exp(logP);
        }) {
    }


    // N.B. Only use this when likelihood is reasonably high
    DOMAIN operator()() {
        DOMAIN sample=priorSampler();
        while(Random::nextDouble() >= likelihood(sample)) {
            sample = priorSampler();
        }
        return sample;
    }

};


#endif //GLPKTEST_REJECTIONSAMPLER_H
