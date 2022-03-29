//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H

#include <functional>
#include "ModelState.h"

template<class DOMAIN>
class RejectionSampler {
public:
    std::function<DOMAIN()>                 priorSampler;
    std::function<double(const DOMAIN &)>   likelihood;

    RejectionSampler(
            std::function<DOMAIN()>                 PriorSampler,
            std::function<double(const DOMAIN &)>   Likelihood
    ) :
     priorSampler(std::move(PriorSampler)),
     likelihood(std::move(Likelihood)) { }


    // N.B. Only use this when likelihood is reasonably high
    DOMAIN operator()() {
        DOMAIN sample=priorSampler();
        while(Random::nextDouble() >= likelihood(sample)) {
//            std::cout << "Rejecting sample " << sample << " with likelihood " << likelihood(sample) << std::endl;
            sample = priorSampler();
        }
//        std::cout << "Accepted sample " << sample << " with likelihood " << likelihood(sample) << std::endl;
        return sample;
    }

};


#endif //GLPKTEST_REJECTIONSAMPLER_H
