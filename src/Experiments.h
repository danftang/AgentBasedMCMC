//
// Created by daniel on 06/05/2021.
//

#ifndef GLPKTEST_EXPERIMENTS_H
#define GLPKTEST_EXPERIMENTS_H

#include <functional>
#include "Trajectory.h"
#include "RejectionSampler.h"
#include "Forecast.h"
#include "Likelihood.h"
#include "diagnostics/AgentDataflow.h"

using namespace dataflow;

class Experiments {
public:
    static void BinomialAgentSingleObservation();
    static void CatMouseSingleObservation();
    static void PredPreySingleObservation();
    static void CatMouseAssimilation();

    template<class AGENT>
    static void doSingleObservationExperiment(
            int nTimesteps,
            int nBurnin,
            int nSamples,
            int nRejectionSamples,
            double kappa,
            const StartStateDistribution<AGENT> &startState,
            const AgentStateObservation<AGENT> &observation);
};


#endif //GLPKTEST_EXPERIMENTS_H
