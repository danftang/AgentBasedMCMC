//
// Created by daniel on 19/08/2021.
//

#ifndef GLPKTEST_MCMCSOLVER_H
#define GLPKTEST_MCMCSOLVER_H

#include "SimplexMCMC.h"
#include "AssimilationWindow.h"
#include "IntSampleStatistics.h"

//template<typename AGENT>
//class MCMCSolver {
//public:
//
//    SimplexMCMC      sampler;
//    IntSampleStatistics exactEndState;
//
//    MCMCSolver(const ConvexPMF &pmf, const std::vector<double> &initialState)
//    :
//    sampler(pmf, initialState),
//    exactEndState(AGENT::domainSize()) {
//    }
//
//    MCMCSolver(const AssimilationWindow<AGENT> &window): MCMCSolver(window.posterior, (window.priorSampler)()) { }
//
//
//    void solve(int nBurnInSamples, int nSamples) {
//        std::cout << "Starting burn-in" << std::endl;
//        for(int burnIn=0; burnIn<nBurnInSamples; ++burnIn) {
//            sampler.nextSample();
//        }
//        std::cout << "Done" << std::endl;
//        //        std::cout << "Sampler:\n" << sampler << std::endl;
//        std::cout << "Initial exactEndState: " << sampler.X() << std::endl;
////        std::cout << "Initial end state: " << Trajectory<AGENT>(sampler.X()).endState() << std::endl;
//        for(int s=0; s<nSamples; ++s) {
//            Trajectory<AGENT> sample(sampler.nextSample());
//            //            std::cout << "Sampler:\n" << sampler << std::endl;
////                        std::cout << "Sample end state: " << sample.endState() << std::endl;
//            exactEndState += sample.endState();
//        }
//        std::cout << "Feasible stats:\n" << sampler.feasibleStatistics << std::endl;
//        std::cout << "Infeasible stats:\n" << sampler.infeasibleStatistics << std::endl;
//        std::cout << "Infeasible proportion = " << sampler.infeasibleStatistics.nSamples*1.0/(sampler.feasibleStatistics.nSamples + sampler.infeasibleStatistics.nSamples) << std::endl;
//    }
//
//    const std::vector<double> &nextSample() {
//        return sampler.nextSample();
//    }
//
//
//
//};


#endif //GLPKTEST_MCMCSOLVER_H
