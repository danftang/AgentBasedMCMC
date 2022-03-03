//
// Created by daniel on 26/08/2021.
//

#include <cmath>
#include <assert.h>
#include "ABMImportanceSampler.h"

std::pair<std::vector<double>, double> ABMImportanceSampler::nextSample() {
    double logWeight;
    std::vector<double> sample;
    do {
        sample = sampler();
        logWeight = logP(sample);
    } while(logWeight == -INFINITY);
    return std::pair(sample, logWeight);
}

std::vector<double> ABMImportanceSampler::expectation(int nSamples) {
    assert(nSamples > 0);
    auto [firstSample, firstLogWeight] = nextSample();
    double logNormalisation = firstLogWeight;
    std::vector<double> mean = std::move(firstSample);
    for(int s=1; s<nSamples; ++s) {
        auto [sample, logWeight] = nextSample();
        assert(sample.size() == mean.size());
        for(int i=0; i<sample.size(); ++i) mean[i] += sample[i]*exp(logWeight - logNormalisation);
    }
    for(int i=0; i<mean.size(); ++i) mean[i] /= exp(logNormalisation);
    return mean;
}
