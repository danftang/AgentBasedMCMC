//
// Created by daniel on 21/10/2021.
//

#ifndef GLPKTEST_VALARRAYLOGGER_H
#define GLPKTEST_VALARRAYLOGGER_H

#include <valarray>
#include <iostream>

template<typename SAMPLE>
class ValarrayLogger {
public:
    std::valarray<SAMPLE> samples;
    int nextFreeIndex;

    ValarrayLogger(int nSamples): samples(nSamples), nextFreeIndex(0) { }

    void operator()(const SAMPLE &sample) {
        if(nextFreeIndex < samples.size()) {
            samples[nextFreeIndex++] = sample;
        } else {
            std::cerr << "Warning: ValarrayLogger full" << std::endl;
        }
    }

    void operator()(SAMPLE &&sample) {
        if(nextFreeIndex < samples.size()) {
            samples[nextFreeIndex++] = std::move(sample);
        } else {
            std::cerr << "Warning: ValarrayLogger full" << std::endl;
        }
    }
};


#endif //GLPKTEST_VALARRAYLOGGER_H
