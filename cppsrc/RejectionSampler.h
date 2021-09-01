//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H


template<typename DOMAIN>
class RejectionSampler {
public:
    const std::function<DOMAIN()> &   priorSampler;
    const ConvexPMF<DOMAIN> &         likelihood;

    RejectionSampler(const std::function<DOMAIN()> &priorSampler,const ConvexPMF<DOMAIN> &likelihood)
    :
    priorSampler(priorSampler),
    likelihood(likelihood) { }

    template<typename AGENT>
    RejectionSampler(const AssimilationWindow<AGENT> &window)
    :
    priorSampler(window.priorSampler),
    likelihood(window.likelihoodPMF)
    { }


    // N.B. Only use this when likelihood is reasonably high
    DOMAIN operator()() {
        DOMAIN sample = priorSampler();
        while(Random::nextDouble() > likelihood.P(sample)) {
//            std::cout << likelihood.P(sample) << " " << sample << std::endl;
            sample = priorSampler();
        }
//        std::cout << "Got sample " << sample << std::endl;
        return sample;
    }

};


#endif //GLPKTEST_REJECTIONSAMPLER_H
