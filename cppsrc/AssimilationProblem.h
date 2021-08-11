//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_ASSIMILATIONPROBLEM_H
#define GLPKTEST_ASSIMILATIONPROBLEM_H

#include "ConvexPMF.h"
#include "ConvexPMFProduct.h"
#include "TrajectoryDistribution.h"
#include "Distribution.h"

class AssimilationProblem: Distribution {
public:
    ConvexPMF                               priorPMF;
    std::function<std::vector<double>()>    priorSampler;
    ConvexPMFProduct                        likelihoodPMF;

    AssimilationProblem(const Distribution &trajectoryPrior)
    : priorPMF(trajectoryPrior.PMF()),
    priorSampler(trajectoryPrior.sampler()),
    likelihoodPMF(priorPMF.nDimensions)
    { }

    AssimilationProblem(ConvexPMF priorPMF, std::function<std::vector<double>()> priorSampler)
    : priorPMF(std::move(priorPMF)),
    priorSampler(std::move(priorSampler)),
    likelihoodPMF(priorPMF.nDimensions)
    {
    }


    void addObservation(ConvexPMF observationLikelihood) {
        likelihoodPMF *= std::move(observationLikelihood);
    }


    SimplexMCMC simplexSampler() const {
        if(likelihoodPMF.size() == 0) return SimplexMCMC(likelihoodPMF * priorPMF, priorSampler());
        return SimplexMCMC(likelihoodPMF * priorPMF, priorSampler());
    }

    Distribution::Sampler sampler() const {
        if(likelihoodPMF.size() == 0) return priorSampler;
        return [mcmc = simplexSampler()]() mutable { return std::vector<double>(mcmc.nextSample()); };
    }

    ConvexPMF PMF() const {
        if(likelihoodPMF.size() == 0) return priorPMF;
        return likelihoodPMF * priorPMF;
    }

    int dimension() const { return priorPMF.nDimensions; }

};


#endif //GLPKTEST_ASSIMILATIONPROBLEM_H
