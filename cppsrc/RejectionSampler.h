//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_REJECTIONSAMPLER_H
#define GLPKTEST_REJECTIONSAMPLER_H

#include "AssimilationProblem.h"

class RejectionSampler {
public:
    const AssimilationProblem &problem;
    RejectionSampler(const AssimilationProblem &prob): problem(prob) {}

    // N.B. Only use this when likelihood is reasonably high
    std::vector<double> operator()() {
        double logLikelihood;
        std::vector<double> sample(0);
        do {
            sample = problem.priorSampler();
            logLikelihood = problem.likelihoodPMF.logP(sample);
        } while(Random::nextDouble() > exp(logLikelihood));
        return sample;
    }

};


#endif //GLPKTEST_REJECTIONSAMPLER_H
