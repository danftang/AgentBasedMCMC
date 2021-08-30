//
// Created by daniel on 25/08/2021.
//

#ifndef GLPKTEST_INTSAMPLESTATISTICS_H
#define GLPKTEST_INTSAMPLESTATISTICS_H

#include <vector>
#include <cmath>
#include "ConvexPMF.h"
#include "Distribution.h"

// Class to represent the marginalised histograms of samples from
// multidimensional integer grids
class IntSampleStatistics: public Distribution {
public:
    std::vector<std::vector<int>> histograms;
    int nSamples;

    IntSampleStatistics(int dimensions)
    : histograms(dimensions),
    nSamples(0) { }

    // pdf takes (dimensionId, count) and returns probability
    IntSampleStatistics(int dimensions, const std::function<double(int, int)> &pdf);

    IntSampleStatistics &operator +=(const std::vector<double> &sample);
    ConvexPMF PMF() const;
    std::function<std::vector<double>()> sampler() const;
    double logP(const std::vector<double> &X) const;
    ConvexPolyhedron constraints() const;
    std::vector<double> nextSample() const;
    int nDimensions() const { return histograms.size(); }

    std::vector<double> means() const;

    double P(int dimension, int count) const {
        if(histograms[dimension].size() > count) return histograms[dimension][count]*1.0/nSamples;
        return 0.0;
    }

    template<typename AGENT>
    static IntSampleStatistics endStateExpectation(const std::function<std::vector<double>()> &trajectorySampler, int nSamples);

    friend std::ostream &operator <<(std::ostream &out, const IntSampleStatistics &sampleStats);
};


#endif //GLPKTEST_INTSAMPLESTATISTICS_H
