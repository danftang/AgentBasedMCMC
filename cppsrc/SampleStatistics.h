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
    int nSamples;

    SampleStatistics(int nDimensions)
    : nSamples(0),
    sum(nDimensions, 0.0),
    sumOfSquares(nDimensions, 0.0)
    { }

    SampleStatistics(const std::function<std::vector<double>()> &sampler, int nSamples): nSamples(nSamples) {
        assert(nSamples > 0);
        std::vector<double> firstSample = sampler();
        sum.resize(firstSample.size(), 0.0);
        sumOfSquares.resize(firstSample.size(), 0.0);
        (*this) += firstSample;
        for(int n=1; n<nSamples; ++n) {
            (*this) += sampler();
        }
    }

    SampleStatistics & operator +=(const std::vector<double> &sample) {
        assert(sample.size() == sum.size());
        ++nSamples;
        for(int i=0; i<nDimensions(); ++i) {
            sum[i] += sample[i];
            sumOfSquares[i] += sample[i]*sample[i];
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

    friend std::ostream &operator <<(std::ostream &out, const SampleStatistics &stats) {
        for(int i=0; i<stats.nDimensions(); ++i) {
            out << "(" << stats.mean(i) << ", " << stats.variance(i) << ") ";
        }
        out << std::endl;
        return out;
    }

};


#endif //GLPKTEST_SAMPLESTATISTICS_H
