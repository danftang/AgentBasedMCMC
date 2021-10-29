//
// Created by daniel on 30/08/2021.
//

#ifndef GLPKTEST_MCMCSAMPLER_H
#define GLPKTEST_MCMCSAMPLER_H

#include "SimplexMCMC.h"

template<typename DOMAIN>
class MCMCSampler {
public:
    SimplexMCMC simplex;
    std::thread *thread;

    MCMCSampler(const ConvexPMF<DOMAIN> &pmf, const DOMAIN &initialState = std::vector<double>())
    : simplex(
        pmf.convexSupport.toLPProblem(),
        [logP = pmf.extendedLogProb](const std::vector<double> &X) { return logP(reinterpret_cast<const DOMAIN &>(X)); },
        initialState
    ),
    thread(NULL) { }

    MCMCSampler(MCMCSampler &&other): simplex(std::move(other.simplex)), thread(other.thread) {
        other.thread = NULL;
        std::cout << "Moving MCMCSampler" << std::endl;
    }

    MCMCSampler(const MCMCSampler<DOMAIN> &other) = delete; // copying a sampler is a bad idea, construct again from pmf


    ~MCMCSampler() {
        if(thread != NULL) delete(thread);
    }

    DOMAIN nextSample() const { return DOMAIN(const_cast<SimplexMCMC &>(simplex).nextSample()); }

    const DOMAIN &operator()() const {
        return static_cast<const DOMAIN &>(const_cast<SimplexMCMC &>(simplex).nextSample());
    }

    // Loggers can be any number of objects that have operator() to consume samples
//    template<typename...LOGGERS>
//    void sample(int nSamples, LOGGERS &... loggers) {
//        for(int s=0; s < nSamples; ++s) {
//            DOMAIN sample = reinterpret_cast<const DOMAIN &>(simplex.nextSample());
//            (loggers(sample),...);
//        }
//    }

    void burnIn(int nSamples) {
        for(int s=0; s<nSamples; ++s) simplex.nextSample();
    }

    void sample(int nSamples, std::function<void(const DOMAIN &)> logger) {
        for(int s=0; s < nSamples; ++s) {
            DOMAIN sample = reinterpret_cast<const DOMAIN &>(simplex.nextSample());
            logger(sample);
        }
    }


//    template<typename...LOGGERS>
//    void sampleInThread(int nSamples, LOGGERS &... loggers) {
//        if(thread == NULL) {
//            thread = new std::thread(&MCMCSampler<DOMAIN>::sample<LOGGERS...>, std::ref(*this), nSamples, loggers...);
////            thread = new std::thread([this, nSamples, &loggers...]() { sample(nSamples, loggers...); });
//        } else {
//            throw("Already running a thread");
//        }
//    }

    void sampleInThread(int nSamples, std::function<void(const DOMAIN &)> logger) {
        if(thread == NULL)
            thread = new std::thread(&MCMCSampler<DOMAIN>::sample, this, nSamples, logger);
         else
             throw("Already running a thread");
    }


    void join() {
        if(thread != NULL) {
            thread->join();
            delete(thread);
            thread = NULL;
        }
    }

    operator std::function<DOMAIN()>() const {
//        std::cout << "Converting MCMCSampler to function" << std::endl;
        return [this]() { return this->nextSample(); };
    }
};


#endif //GLPKTEST_MCMCSAMPLER_H
