//
// Created by daniel on 18/08/2021.
//

#ifndef GLPKTEST_TRAJECTORYSAMPLER_H
#define GLPKTEST_TRAJECTORYSAMPLER_H

template<typename AGENT>
class TrajectorySampler {
public:
    int nTimesteps;
    std::function<Trajectory<AGENT>()> sampler;

    TrajectorySampler(int nTimesteps, const std::function<std::vector<double>()> &startStateSampler)
    : nTimesteps(nTimesteps),
    sampler([startStateSampler,nTimesteps]() {
        return Trajectory<AGENT>(nTimesteps, startStateSampler);
    }) { }

    std::vector<double> operator()() const {
        return sampler();
    }


    int nDimensions() { return Trajectory<AGENT>::dimension(nTimesteps); }


    Trajectory<AGENT> nextSample() const {
        return sampler();
    }

    std::function<std::vector<double>()> endStateSampler() {
        return [*this]() {
            return nextSample()(nTimesteps);
        };
    }

    IntSampleStatistics endState(int nSamples) {
        IntSampleStatistics stats(AGENT::domainSize());
        for(int s=0; s<nSamples; ++s) {
//            std::cout << "prior end state = " << nextSample()(nTimesteps);
//            stats += nextSample().endState();//(nTimesteps);
            nextSample();//.endState();//(nTimesteps);
        }
        return stats;
    }
};


#endif //GLPKTEST_TRAJECTORYSAMPLER_H
