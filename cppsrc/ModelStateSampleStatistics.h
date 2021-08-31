//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_MODELSTATESAMPLESTATISTICS_H
#define GLPKTEST_MODELSTATESAMPLESTATISTICS_H

#include "IntSampleStatistics.h"

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
        for(int s = 0; s<nSamples; ++s) {
            ModelState<AGENT> endState = trajectorySampler().endState();
//            std::cout << "Taking sample " << s << " " << endState << std::endl;
            (*this) += endState;
        }
//        std::cout << "Done sampling with " << this->means() << std::endl;
    }

};


#endif //GLPKTEST_MODELSTATESAMPLESTATISTICS_H
