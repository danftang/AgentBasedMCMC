//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H


class RejectionSampler {
public:
    const std::function<std::vector<double>()> &   priorSampler;
    const ConvexPMF &                              likelihood;

    RejectionSampler(const std::function<std::vector<double>()> &priorSampler,const ConvexPMF &likelihood)
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
    std::vector<double> operator()() {
        std::vector<double> sample(0);
        do {
            sample = priorSampler();
        } while(likelihood.convexSupport.isValidSolution(sample) &&  Random::nextDouble() > exp(likelihood.logP(sample)));
        return sample;
    }

};


#endif //GLPKTEST_REJECTIONSAMPLER_H
