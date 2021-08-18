//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H


class RejectionSampler {
public:
    std::function<std::vector<double>()>    priorSampler;
    ConvexPMF                               likelihoodPMF;

    RejectionSampler(std::function<std::vector<double>()> priorSampler, ConvexPMF likelihood)
    : priorSampler(std::move(priorSampler)),
      likelihoodPMF(std::move(likelihood))
      {}

    // N.B. Only use this when likelihood is reasonably high
    std::vector<double> operator()() {
        double logLikelihood;
        std::vector<double> sample(0);
        do {
            sample = priorSampler();
            logLikelihood = likelihoodPMF.logP(sample);
        } while(Random::nextDouble() > exp(logLikelihood));
        return sample;
    }

};


#endif //GLPKTEST_REJECTIONSAMPLER_H
