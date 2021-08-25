//
// Created by daniel on 25/08/2021.
//

#ifndef GLPKTEST_INTSAMPLESTATISTICS_H
#define GLPKTEST_INTSAMPLESTATISTICS_H

#include <vector>
#include <cmath>
#include "ConvexPMF.h"

// Class to represent the marginalised histograms of samples from
// multidimensional integer grids
class IntSampleStatistics {
    std::vector<std::vector<int>> histograms;
    int nSamples;

    IntSampleStatistics &operator +=(const std::vector<double> &sample);
    ConvexPMF PMF() const;
    double logP(const std::vector<double> &X) const;
    ConvexPolyhedron constraints() const;
    int nDimensions() const { return histograms.size(); }
};


#endif //GLPKTEST_INTSAMPLESTATISTICS_H
