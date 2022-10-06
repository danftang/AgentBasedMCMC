//
// Created by daniel on 06/05/2021.
//

#ifndef GLPKTEST_EXPERIMENTS_H
#define GLPKTEST_EXPERIMENTS_H

#include "ConstrainedFactorisedDistribution.h"
#include "ABMPosterior.h"

class Experiments {
public:
    static void BinomialAgentSingleObservation();
    static void CatMouseSingleObservation();
    static void PredPreySingleObservation();
    static void CatMouseAssimilation();
    static void Animation();

    template<class DOMAIN, class STARTSTATE>
    static void doValidationExperiment(
            const ABMPosterior<DOMAIN,STARTSTATE> &,
            int nBurnin,
            int nSamples,
            int nRejectionSamples);
};


#endif //GLPKTEST_EXPERIMENTS_H
