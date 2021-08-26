//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_ABMIMPORTANCESAMPLER_H
#define GLPKTEST_ABMIMPORTANCESAMPLER_H

#include <vector>
#include <functional>

// PROBLEM should be an Assimilation Problem with traits and members as in AssimilationProblem
class ABMImportanceSampler {
public:
    std::function<std::vector<double>()>        sampler;
    std::function<double(const std::vector<double> &)>  logP;

    ABMImportanceSampler(std::function<std::vector<double>()> sampler, std::function<double(const std::vector<double> &)> logP)
    : sampler(sampler),
    logP(logP) {}

    // N.B. Weight is log and unnormalised
    std::pair<std::vector<double>,double> nextSample();

    std::vector<double> expectation(int nSamples);
};


#endif //GLPKTEST_ABMIMPORTANCESAMPLER_H
