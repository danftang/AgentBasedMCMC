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
#include "../FactorisedDistributionOnBasis.h"
#include "../FactorisedDistributionSampler.h"
#include "../ModelState.h"

class ChainStats {
public:
    MeanAndVariance meanVariance;
    std::valarray<std::valarray<double>> variogram;
    int varioStride;
    std::valarray<double> meanEndState;
    MCMCStatistics samplerStats;

    ChainStats() {}

    template<typename SYNOPSIS>
    ChainStats(const std::valarray<SYNOPSIS> &synopsisSamples,
               int nLags,
               double maxLagProportion,
               std::valarray<double> meanEndState,
               const MCMCStatistics &samplerStats) :
            meanVariance(synopsisSamples),
            variogram(calcVariogram(synopsisSamples, nLags, maxLagProportion)),
            varioStride((maxLagProportion * synopsisSamples.size()) / (nLags - 1.0)),
            meanEndState(std::move(meanEndState)),
            samplerStats(samplerStats)
            { }

    int nSamples() const { return meanVariance.nSamples; }

    int dataDimension() const { return meanVariance.dataDimension(); }

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
    calcVariogram(const std::valarray<SAMPLE> &samples, int nLags = 100, double maxLagProportion = 0.5) {
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

    // Following
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // we want to compare the statistics of the first half and last half of a Markov chain in order
    // to better detect convergence, so this takes samples, and generates first/last half statistics.
    template<class SAMPLER>
    static std::pair<ChainStats,ChainStats> sampleChainStatsPair(SAMPLER &sampler, int nSamples) {
        std::cout << "Starting ChainStats on sampler " << sampler.basisVectors.size() << " " << sampler.factors.size() << std::endl;
        assert(nSamples%1 == 0); // nSamples should be an even number
        const double maxLagProportion = 0.5;
        const int nLags = 200;
        const int nBurnIn = nSamples/5;
        const int endStateDimension = sampler().endState().size();

        std::valarray<std::valarray<double>> firstHalfSynopsisSamples(nSamples / 2);
        std::valarray<std::valarray<double>> lastHalfSynopsisSamples(nSamples / 2);
        std::valarray<double> firstHalfMeanEndState(0.0, endStateDimension);
        std::valarray<double> lastHalfMeanEndState(0.0, endStateDimension);

        for(int s = 1; s<nBurnIn; ++s) sampler();
        for (int s = 0; s < nSamples; ++s) {
            auto endState = sampler().endState();
            if (s < nSamples / 2) {
                firstHalfSynopsisSamples[s] = convergenceStatistics(endState);
                firstHalfMeanEndState += endState;
            } else {
                lastHalfSynopsisSamples[s - nSamples/2] = convergenceStatistics(endState);
                lastHalfMeanEndState += endState;
            }
        }

        std::cout << "Stats =\n" << sampler.stats << std::endl;

        firstHalfMeanEndState /= nSamples/2.0;
        lastHalfMeanEndState /= nSamples/2.0;

        return std::pair(
                ChainStats(firstHalfSynopsisSamples, nLags, maxLagProportion, firstHalfMeanEndState, sampler.stats),
                ChainStats(lastHalfSynopsisSamples, nLags, maxLagProportion, lastHalfMeanEndState, sampler.stats));
    }


private:
    friend class boost::serialization::access;

    template<typename T> static bool hasGridsize() { return hasGridsize<T>(0); }
    template<typename T, int = T::gridsize> static bool hasGridsize(int dummy) { return true; }
    template<typename T, typename DUMMYARG> static bool hasGridsize(DUMMYARG x) { return false; }

    // For gridded agents, use total occupations of decreasing size square areas
    // aligned along the diagonal
    template<class AGENT, int GRIDSIZE = AGENT::gridsize>
    static std::valarray<double> convergenceStatistics(const ModelState<AGENT> &endState) {
        std::valarray<double> synopsis(floor(log2(GRIDSIZE)) -1);
        int origin = 0;
        int varid = 0;
        for(int partitionSize = GRIDSIZE / 2; partitionSize > 1; partitionSize /=2) {
            int occupation = 0;
            for(int x=0; x < partitionSize; ++x) {
                for(int y=0; y < partitionSize; ++y) {
                    for(int type = 0; type <= AGENT::typeDomainSize; ++type) {
                        occupation += endState[AGENT(origin + x,origin + y,typename AGENT::Type(type))];
                    }
                }
            }
            assert(varid < synopsis.size());
            synopsis[varid++] = occupation;
            origin += partitionSize;
        }
        return synopsis;
    }

    template<class AGENT, typename = std::enable_if_t<!hasGridsize<AGENT>()>>
    static std::valarray<double> convergenceStatistics(const ModelState<AGENT> &endState) {
        std::valarray<double> synopsis(floor(log2(AGENT::domainSize)) -1);
        int origin = 0;
        int varid = 0;
        for(int partitionSize = AGENT::domainSize / 2; partitionSize > 1; partitionSize /=2) {
            int occupation = 0;
            for(int agentId = origin; agentId < origin + partitionSize; ++agentId) {
                occupation += endState[AGENT(agentId)];
            }
            assert(varid < synopsis.size());
            synopsis[varid++] = occupation;
            origin += partitionSize;
        }
        return synopsis;
    }




    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & meanVariance & variogram & varioStride & samplerStats;
    }
};



#endif //GLPKTEST_CHAINSTATS_H
