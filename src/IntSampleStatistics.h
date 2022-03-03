//
// Created by daniel on 25/08/2021.
//

#ifndef GLPKTEST_INTSAMPLESTATISTICS_H
#define GLPKTEST_INTSAMPLESTATISTICS_H

#include <vector>
#include <cmath>
#include "ConvexPMF.h"
#include "Distribution.h"
#include "StlStream.h"

// Class to represent the marginalised histograms of samples from
// multidimensional integer grids
template<typename DOMAIN>
class IntSampleStatistics: public Distribution<DOMAIN> {
public:
    std::vector<std::vector<int>>   histograms;
    int nSamples;

    IntSampleStatistics(int dimensions)
    : histograms(dimensions),
    nSamples(0) { }

    int nDimensions() const { return histograms.size(); }

    // pdf takes (dimensionId, count) and returns probability
    IntSampleStatistics(int dimensions, const std::function<double(int, int)> &pdf): IntSampleStatistics(dimensions) {
        nSamples = 1000000;
        for(int i=0; i<nDimensions(); ++i) {
            double cumulativeP = 0.0;
            int cumulativeCount = 0;
            int j = 0;
            while(cumulativeCount < nSamples) {
                cumulativeP += pdf(i,j);
                int count = std::round(cumulativeP*nSamples) - cumulativeCount;
                cumulativeCount += count;
                if(cumulativeCount > nSamples) count -= cumulativeCount - nSamples;
                histograms[i].push_back(count);
                ++j;
            }
        }
    }

    IntSampleStatistics &operator +=(const DOMAIN &sample) {
        ++nSamples;
        if(histograms.size() < sample.size()) histograms.resize(sample.size());
        for(int i=0; i<sample.size(); ++i) {
            int Xi = std::round(sample[i]);
            if(histograms[i].size() < Xi + 1) histograms[i].resize(Xi + 1, 0);
            ++(histograms[i][Xi]);
        }
        return *this;
    }

//    void sampleFrom(const std::function<DOMAIN()> &sampler, int nSamples) {
//        for(int s=0; s<nSamples; ++s) {
//            (*this) += sampler();
//        }
//    }


    ConvexPMF<DOMAIN> PMF() const {
        return ConvexPMF<DOMAIN>(
            [*this, lBounds = lowerBounds(), logMargs = logMarginals()](const DOMAIN &X) {
                double logP = extendedLogP(X,lBounds,logMargs);
                return logP;
            },
            nDimensions(),
            constraints()
        );
    }

    std::function<DOMAIN()> sampler() const { return *this; }


    double logP(const DOMAIN &X) const {
        double lp = 0.0;
        for(int i=0; i<X.size(); ++i) {
            int Xi = std::round(X[i]);
            const std::vector<int> &histogram = histograms[i];
            if(Xi < histogram.size()) {
                lp += log(histogram[Xi]*1.0/nSamples);
            } else {
                lp -= INFINITY;
                i = X.size();
            }
        }
        return lp;
    }

    double extendedLogP(const DOMAIN &X, const std::vector<int> &lowerBounds, const std::vector<double> &logMarginals) const {
        double lp = 0.0;
        for(int i=0; i<X.size(); ++i) {
            int Xi = std::round(X[i]);
            const std::vector<int> &histogram = histograms[i];
            if(lowerBounds[i] <= Xi && Xi < histogram.size()) {
                lp += log(histogram[Xi]*1.0/nSamples);
            } else {
                lp += logMarginals[i];
            }
        }
        return lp;
    }


    ConvexPolyhedron constraints() const {
        ConvexPolyhedron constraints;
        assert(nSamples > 0);
        constraints.reserve(nDimensions());
        std::vector<int> lBounds = lowerBounds();
        for(int i=0; i<nDimensions(); ++i) {
            constraints.push_back(lBounds[i] <= 1.0*glp::X(i) <= histograms.size()-1);
        }
        return constraints;
    }


    std::vector<int> lowerBounds() const {
        std::vector<int> lBounds(nDimensions());
        assert(nSamples > 0);
        for(int i=0; i<nDimensions(); ++i) {
            int lowerBound = 0;
            while(histograms[i][lowerBound] == 0) {
                ++lowerBound;
                assert(lowerBound < histograms[i].size());
            }
            lBounds[i] = lowerBound;
        }
        return lBounds;
    }


    std::vector<double> logMarginals() const {
        assert(nSamples > 0);
        std::vector<double> logMargs(histograms.size());
        for(int i=0; i<histograms.size(); ++i) {
            double expectationP = 0.0;
            const std::vector<int> &histogram = histograms[i];
            for (int j = 0; j < histogram.size(); ++j) {
                double p = histogram[j]*1.0/nSamples;
                expectationP += p*p;
            }
            logMargs[i] = log(infeasibleExpectationFraction * expectationP);
        }
//        std::cout << "Log marginals " << logMargs << std::endl;
        return logMargs;
    }


    DOMAIN operator()() const { return nextSample(); }
    DOMAIN nextSample() const {
        std::vector<double> sample(nDimensions());
        for(int i=0; i<nDimensions(); ++i) {
            int targetCumulativeCount = Random::nextInt(0, nSamples);
            int cumulativeCount = 0;
            int j=-1;
            do {
                cumulativeCount += histograms[i][++j];
                assert(j < histograms[i].size());
            } while(cumulativeCount <= targetCumulativeCount);
            sample[i] = j;
        }
        return sample;
    }


    void clear() {
        for(int i=0; i<histograms.size(); ++i) histograms[i].clear();
        nSamples = 0;
    }


    DOMAIN means() const {
        std::vector<double> mu(nDimensions());
        for(int i=0; i<nDimensions(); ++i) {
            mu[i] = 0.0;
            for(int j=1; j<histograms[i].size(); ++j) {
                mu[i] += j*histograms[i][j];
            }
            mu[i] /= nSamples;
        }
        return mu;

    }

    double P(int dimension, int count) const {
        if(histograms[dimension].size() > count) return histograms[dimension][count]*1.0/nSamples;
        return 0.0;
    }

//    template<typename AGENT>
//    static IntSampleStatistics endStateExpectation(const std::function<std::vector<double>()> &trajectorySampler, int nSamples);

    friend std::ostream &operator <<(std::ostream &out, const IntSampleStatistics &sampleStats) {
        for(const std::vector<int> &histogram: sampleStats.histograms) {
            out << histogram << std::endl;
        }
        return out;
    }
};


#endif //GLPKTEST_INTSAMPLESTATISTICS_H
