//
// Created by daniel on 25/08/2021.
//

#include "IntSampleStatistics.h"
#include "Random.h"

IntSampleStatistics::IntSampleStatistics(int dimensions, const std::function<double(int, int)> &pdf) : IntSampleStatistics(dimensions) {
    nSamples = 1000000;
    for(int i=0; i<nDimensions(); ++i) {
        double cumulativeP = 0.0;
        int cumulativeCount = 0;
        int j = 0;
        while(cumulativeCount < nSamples) {
            cumulativeP += pdf(i,j);;
            int count = std::round(cumulativeP*nSamples) - cumulativeCount;
            cumulativeCount += count;
            if(cumulativeCount > nSamples) count -= cumulativeCount - nSamples;
            histograms[i].push_back(count);
            ++j;
        }
    }
}

IntSampleStatistics &IntSampleStatistics::operator+=(const std::vector<double> &sample) {
    ++nSamples;
    if(histograms.size() < sample.size()) histograms.resize(sample.size());
    for(int i=0; i<sample.size(); ++i) {
        int Xi = std::round(sample[i]);
        if(histograms[i].size() < Xi + 1) histograms[i].resize(Xi + 1, 0);
        ++(histograms[i][Xi]);
    }
    return *this;
}

ConvexPMF IntSampleStatistics::PMF() const {
    return ConvexPMF([*this](const std::vector<double> &X) {
                         return logP(X);
                     },
                     nDimensions(),
                     constraints());
}

double IntSampleStatistics::logP(const std::vector<double> &X) const {
    double lp = 0.0;
    for(int i=0; i<X.size(); ++i) {
        int Xi = std::round(X[i]);
        const std::vector<int> &histogram = histograms[i];
        if(histogram.size() > Xi) {
            lp += log(histogram[Xi]*1.0/nSamples);
        } else {
            lp -= INFINITY;
            i = X.size();
        }
    }
    return lp;
}

ConvexPolyhedron IntSampleStatistics::constraints() const {
    ConvexPolyhedron constraints;
    assert(nSamples > 0);
    for(int i=0; i<nDimensions(); ++i) {
        int lowerBound = 0;
        while(histograms[i][lowerBound] == 0) {
            ++lowerBound;
            assert(lowerBound < histograms[i].size());
        }
        constraints.push_back(lowerBound <= 1.0*glp::X(i) <= histograms.size()-1);
    }
    return constraints;
}

std::function<std::vector<double>()> IntSampleStatistics::sampler() const {
    return [*this]() { return nextSample(); };
}

std::vector<double> IntSampleStatistics::nextSample() const {
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

std::vector<double> IntSampleStatistics::means() const {
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


std::ostream &operator <<(std::ostream &out, const IntSampleStatistics &sampleStats) {
    for(const std::vector<int> &histogram: sampleStats.histograms) {
        out << histogram << std::endl;
    }
    return out;
}
