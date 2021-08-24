//
// Created by daniel on 18/08/2021.
//

#ifndef GLPKTEST_TRAJECTORYPRIORSAMPLER_H
#define GLPKTEST_TRAJECTORYPRIORSAMPLER_H

template<typename AGENT>
class TrajectoryPriorSampler {
public:
    int nTimesteps;
    std::function<std::vector<double>()> startStateSampler;

    TrajectoryPriorSampler(int nTimesteps, std::function<std::vector<double>()> startStateSampler)
    : nTimesteps(nTimesteps),
    startStateSampler(startStateSampler) {}

    std::vector<double> operator()() const {
        return Trajectory<AGENT>(nTimesteps, startStateSampler);
    }

    std::function<std::vector<double>()> endStateSampler() {
        return [*this]() {
            return Trajectory<AGENT>(nTimesteps, startStateSampler)(nTimesteps);
        };
    }

    SampleStatistics endState(int nSamples) {
        SampleStatistics stats(AGENT::domainSize());
        for(int s=0; s<nSamples; ++s) {
            stats += Trajectory<AGENT>(nTimesteps, startStateSampler)(nTimesteps);
        }
        return stats;
    }
};


#endif //GLPKTEST_TRAJECTORYPRIORSAMPLER_H
