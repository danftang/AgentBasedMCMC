//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H

#include <functional>
#include "ModelState.h"
#include "ABMPrior.h"
#include "Likelihood.h"

template<class DOMAIN>
class RejectionSampler {
public:
    std::function<const DOMAIN &()>         priorSampler;
    std::function<double(const DOMAIN &)>   likelihood;

    RejectionSampler(std::function<DOMAIN()> PriorSampler, std::function<double(const DOMAIN &)> Likelihood):
        priorSampler(std::move(PriorSampler)),
        likelihood(std::move(Likelihood)) {
    }


//    RejectionSampler(ABMPrior<DOMAIN> & prior, const ConstrainedFactorisedDistribution<DOMAIN> & likelihood) :
//        priorSampler(prior.sampler()),
//        likelihood([&likelihood](const DOMAIN &X) {
//            return likelihood.Pexact(X);
//        }) {
//    }


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
