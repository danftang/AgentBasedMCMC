//
// Created by daniel on 03/11/2021.
//

#ifndef GLPKTEST_CHAINSTATS_H
#define GLPKTEST_CHAINSTATS_H

#include <valarray>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/valarray.hpp>
#include <boost/serialization/vector.hpp>

#include "MeanAndVariance.h"
#include "../MCMCStatistics.h"

class ChainStats {
public:
    MeanAndVariance meanVariance;
    std::valarray<std::valarray<double>> vario;
    int varioStride;
    std::valarray<double> meanEndState;
    MCMCStatistics stats;

    ChainStats() {}

    template<typename SYNOPSIS>
    ChainStats(std::valarray<SYNOPSIS> synopsisSamples,
               int nLags,
               double maxLagProportion,
               std::valarray<double> meanEndState,
               const MCMCStatistics &Stats) :
            meanVariance(std::move(synopsisSamples)),
            vario(variogram(synopsisSamples, nLags, maxLagProportion)),
            varioStride((maxLagProportion * synopsisSamples.size()) / (nLags - 1.0)),
            meanEndState(std::move(meanEndState)),
            stats(Stats) {
    }

    int nSamples() const { return meanVariance.nSamples; }

    int dimension() const { return meanVariance.dimension(); }

    ////////////////////////////////////////////////////////////////////////////
    // Calculates the variogram, V_t, of discrete samples x_0...x_n as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, page 286.
    //
    // V_t = 1/(n-t) sum_{i=0}^{n-t-1} (x_i - x_{i+t})^2
    //
    // The lag, t, is calculated at "nLags" different values
    // evently spaced from 0 until "maxLagProportion"*(n+1).
    //
    // Since we're only interested in estimating the overall rate of decay,
    // we don't need to calculate g_t at all values of t, so although this
    // could be done by using FFT on the zero-padded samples to convert into
    // and back from the frequency domain, the saving is negligible for the sizes
    // we're interested in (nLags approx 100) so we just calculate explicitly.
    //
    // The SAMPLE type can be any type that has arithmetic operations, notably
    // another valarray.
    //
    ////////////////////////////////////////////////////////////////////////////
    template<typename SAMPLE>
    static std::valarray<std::valarray<double>>
    variogram(const std::valarray<SAMPLE> &samples, int nLags = 100, double maxLagProportion = 0.5) {
        std::valarray<std::valarray<double>> vario(0.0 * samples[0], nLags);

        int stride = (maxLagProportion * samples.size()) / (nLags - 1.0);
        for (int j = 0; j < nLags; ++j) {
            int t = j * stride;
            for (int i = 0; i < samples.size() - t; ++i) {
                auto dxt = samples[i] - samples[i + t];
                vario[j] += dxt * dxt;
            }
            vario[j] /= samples.size() - t;
        }
        return vario;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & meanVariance & vario & varioStride & meanEndState & stats;
    }
};



#endif //GLPKTEST_CHAINSTATS_H
