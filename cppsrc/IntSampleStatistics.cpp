//
// Created by daniel on 25/08/2021.
//

#include "IntSampleStatistics.h"

IntSampleStatistics &IntSampleStatistics::operator+=(const std::vector<double> &sample) {
    ++nSamples;
    if(histograms.size() < sample.size()) histograms.resize(sample.size());
    for(int i=0; i<sample.size(); ++i) {
        int Xi = std::round(sample[i]);
        if(histograms[i].size() < Xi + 1) histograms[i].resize(Xi + 1, 0);
        ++histograms[i][Xi];
    }
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
            lp += log(histogram[Xi]);
        } else {
            lp -= INFINITY;
            i = X.size();
        }
    }
    lp -= X.size()*log(nSamples);
    return lp;
}

ConvexPolyhedron IntSampleStatistics::constraints() const {
    ConvexPolyhedron constraints;
    for(int i=0; i<nDimensions(); ++i) {
        int lowerBound = 0;
        while(lowerBound < histograms[i].size() && histograms[i][lowerBound] == 0) ++lowerBound;
        constraints.push_back(lowerBound <= 1.0*glp::X(i) <= histograms.size()-1);
    }
    return constraints;
}
