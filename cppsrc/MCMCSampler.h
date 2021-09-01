//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_MCMCSAMPLER_H
#define GLPKTEST_MCMCSAMPLER_H

#include "SimplexMCMC.h"

template<typename DOMAIN>
class MCMCSampler {
    MCMCSampler(const MCMCSampler &other): simplex(other.simplex) {
        std::cout << "Copying MCMCSampler" << std::endl;
    }

public:
    SimplexMCMC simplex;

    MCMCSampler(MCMCSampler &&other): simplex(std::move(other.simplex)) {
        std::cout << "Moving MCMCSampler" << std::endl;
    }

    MCMCSampler(const ConvexPMF<DOMAIN> &pmf, const DOMAIN &initialState = std::vector<double>())
    : simplex(
        pmf.convexSupport.toLPProblem(),
        [logP = pmf.extendedLogProb](const std::vector<double> &X) { return logP(reinterpret_cast<const DOMAIN &>(X)); },
        initialState
    ) { }

    DOMAIN nextSample() const { return DOMAIN(const_cast<SimplexMCMC &>(simplex).nextSample()); }

    operator std::function<DOMAIN()>() const {
        std::cout << "Converting MCMCSampler to function" << std::endl;
        return [this]() { return this->nextSample(); };
    }
};


#endif //GLPKTEST_MCMCSAMPLER_H
