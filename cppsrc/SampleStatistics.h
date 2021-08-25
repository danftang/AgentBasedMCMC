//
// Created by daniel on 17/08/2021.
//

#ifndef GLPKTEST_SAMPLESTATISTICS_H
#define GLPKTEST_SAMPLESTATISTICS_H

#include "valarray"

// Represents statistics of a set of sample vectors
class SampleStatistics {
public:
    std::valarray<double> sum;
    std::valarray<double> sumOfSquares;
    std::valarray<double> max;
    int nSamples;

    SampleStatistics(int nDimensions)
    : nSamples(0),
    sum(0.5, nDimensions),
    sumOfSquares(0.25, nDimensions),
    max(-INFINITY, nDimensions)
    { }

    SampleStatistics(const std::function<std::vector<double>()> &sampler, int nSamples): nSamples(0) {
        assert(nSamples > 0);
        std::vector<double> firstSample = sampler();
        sum.resize(firstSample.size(), 0.5);
        sumOfSquares.resize(firstSample.size(), 0.25);
        max.resize(firstSample.size(), -INFINITY);
        (*this) += firstSample;
        for(int n=1; n<nSamples; ++n) {
            (*this) += sampler();
        }
    }

    SampleStatistics & operator +=(const std::vector<double> &sample) {
        assert(sample.size() == nDimensions());
        ++nSamples;
        for(int i=0; i<nDimensions(); ++i) {
            double s = sample[i];
            sum[i] += s;
            sumOfSquares[i] += s*s;
            if(s > max[i]) max[i] = s;
        }
        return *this;
    }

    double mean(int i) const { return sum[i]/nSamples; }
    double variance(int i) const { return (sumOfSquares[i] - sum[i]*sum[i]/nSamples)/nSamples; }

    std::valarray<double> means() const {
        return sum * (1.0/nSamples);
    }

    std::valarray<double> variances() const {
        return (sumOfSquares - sum*sum*(1.0/ nSamples))*(1.0/nSamples);
    }

    int nDimensions() const { return sum.size(); }

    void reset() {
        nSamples = 0;
        sum = 0.5;
        sumOfSquares = 0.25;
        max = -INFINITY;
    }

    friend std::ostream &operator <<(std::ostream &out, const SampleStatistics &stats) {
        for(int i=0; i<stats.nDimensions(); ++i) {
            out << "(" << stats.mean(i) << ", " << stats.variance(i) << ") ";
        }
        out << std::endl;
        return out;
    }

};


#endif //GLPKTEST_SAMPLESTATISTICS_H
