//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_MODELSTATESAMPLESTATISTICS_H
#define GLPKTEST_MODELSTATESAMPLESTATISTICS_H

#include "IntSampleStatistics.h"
#include "diagnostics/KLDivergence.h"

template<typename AGENT>
class ModelStateSampleStatistics: public IntSampleStatistics<ModelState<AGENT>> {
public:
    ModelStateSampleStatistics()
    : IntSampleStatistics<ModelState<AGENT>>(AGENT::domainSize()) { }

    ModelStateSampleStatistics(const std::function<Trajectory<AGENT>()> &trajectorySampler, int nSamples)
    : ModelStateSampleStatistics() {
        sampleFromEndState(trajectorySampler, nSamples);
    }

    void sampleFromEndState(const std::function<Trajectory<AGENT>()> &trajectorySampler, int nSamples) {
//        this->clear();
        KLDivergence divergenceStats(nSamples/100, nSamples*0.9/200);
        for(int s = 0; s<nSamples; ++s) {
            Trajectory<AGENT> T = trajectorySampler();
            divergenceStats << T.logProb();
            ModelState<AGENT> endState = T.endState();
            if(s%0x10000 == 0) std::cout << "Taking sample " << s << " " << std::endl;//endState << std::endl;
            (*this) += endState;
        }
        divergenceStats.plot();
//        debug();
//        std::cout << "Done sampling with " << this->means() << std::endl;
    }

};


#endif //GLPKTEST_MODELSTATESAMPLESTATISTICS_H
