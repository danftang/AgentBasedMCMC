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
    static constexpr int defaultNLags = 200;
    static constexpr double defaultMaxLagProportion = 0.5;

    MeanAndVariance meanVariance;
    int varioStride;
    std::valarray<std::valarray<double>> variogram;
    std::valarray<double> meanEndState;
//    MCMCStatistics samplerStats;
    std::chrono::steady_clock::duration sampleTime;

    ChainStats() {}

    template<class SAMPLER>
    ChainStats(int nSamples,SAMPLER &sampler, int nLags = defaultNLags, double maxLagProportion = defaultMaxLagProportion):
            varioStride(std::max(1,static_cast<int>(maxLagProportion * nSamples) / (nLags - 1))),
            meanEndState(0.0, sampler().endState().size())
    {

        std::valarray<std::valarray<double>> statisticsSamples(nSamples);

        auto startTime = std::chrono::steady_clock::now();

        for (int s = 0; s < nSamples; ++s) {
            auto endState = sampler().endState();
            statisticsSamples[s] = convergenceStatistics(endState);
            meanEndState += endState;
        }

        auto endTime = std::chrono::steady_clock::now();

//        std::cout << "Stats =\n" << sampler.stats << std::endl;

        meanEndState /= nSamples;
        meanVariance.template addDatapoints(statisticsSamples);
        variogram = calcVariogram(statisticsSamples, nSamples*maxLagProportion, varioStride);
//        samplerStats = sampler.stats;
        sampleTime = endTime - startTime;

    }


//    template<typename SYNOPSIS>
//    ChainStats(const std::valarray<SYNOPSIS> &synopsisSamples,
//               int nLags,
//               int lagStride,
//               std::valarray<double> meanEndState,
//               const MCMCStatistics &samplerStats) :
//            meanVariance(synopsisSamples),
//            varioStride(lagStride),
//            variogram(calcVariogram(synopsisSamples, nLags, varioStride)),
//            meanEndState(std::move(meanEndState)),
//            samplerStats(samplerStats)
//            { }

    int nSamples() const { return meanVariance.nSamples; }

    int dataDimension() const { return meanVariance.dataDimension(); }

    ////////////////////////////////////////////////////////////////////////////
    // Calculates the variogram, V_t, of discrete samples x_0...x_n as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, page 286.
    //
    // V_t = 1/(n-t) sum_{i=0}^{n-t-1} (x_i - x_{i+t})^2
    //
    // The lag, t, is calculated at equally spaced intervals from 0 until (but not including)
    // "maxLag" with a spacing of 'lagStride'
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
    calcVariogram(const std::valarray<SAMPLE> &samples, int maxLag, int lagStride) {
        assert(maxLag <= samples.size());
        std::valarray<std::valarray<double>> vario(0.0 * samples[0], maxLag/lagStride);

        int j = 0;
        for(int lag = 0; lag < maxLag; lag += lagStride) {
            for (int i = 0; i < samples.size() - lag; ++i) {
                auto dxt = samples[i] - samples[i + lag];
                vario[j] += dxt * dxt;
            }
            vario[j] /= samples.size() - lag;
            ++j;
        }
        return vario;
    }

    // Following
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // we want to compare the statistics of the first half and last half of a Markov chain in order
    // to better detect convergence, so this takes samples, and generates first/last half statistics.
//    template<class SAMPLER, class SAMPLE = std::result_of<SAMPLER()>>
//    static std::pair<ChainStats,ChainStats> makeChainStatsPair(SAMPLER &sampler, int nSamples) {
//        assert(nSamples%1 == 0); // nSamples should be an even number
//        const double maxLagProportion = 0.5;
//        const int nLags = 200;
//        const int nBurnIn = nSamples/5;
//        const int endStateDimension = sampler().endState().size();
//
//        std::valarray<std::valarray<double>> firstHalfSynopsisSamples(nSamples / 2);
//        std::valarray<std::valarray<double>> lastHalfSynopsisSamples(nSamples / 2);
//        std::valarray<double> firstHalfMeanEndState(0.0, endStateDimension);
//        std::valarray<double> lastHalfMeanEndState(0.0, endStateDimension);
//
//        for(int s = 1; s<nBurnIn; ++s) sampler();
//        for (int s = 0; s < nSamples; ++s) {
//            auto endState = sampler().endState();
//            if (s < (nSamples / 2)) {
//                firstHalfSynopsisSamples[s] = convergenceStatistics(endState);
//                firstHalfMeanEndState += endState;
//            } else {
//                lastHalfSynopsisSamples[s - nSamples/2] = convergenceStatistics(endState);
//                lastHalfMeanEndState += endState;
//            }
//        }
//
//        std::cout << "Stats =\n" << sampler.stats << std::endl;
//
//        firstHalfMeanEndState /= nSamples/2.0;
//        lastHalfMeanEndState /= nSamples/2.0;
//
////        std::cout << "First half mean end state = " << firstHalfMeanEndState << std::endl;
////        std::cout << "Last  half mean end state = " << lastHalfMeanEndState << std::endl;
//
//        return std::pair(
//                ChainStats(firstHalfSynopsisSamples, nLags, maxLagProportion, firstHalfMeanEndState, sampler.stats),
//                ChainStats(lastHalfSynopsisSamples, nLags, maxLagProportion, lastHalfMeanEndState, sampler.stats));
//    }

    template<class DISTRIBUTION, class TEST = std::enable_if_t<!std::is_invocable_v<DISTRIBUTION>>>
    static std::pair<ChainStats,ChainStats> makeChainStatsPair(const DISTRIBUTION &distribution, int nSamples) {
        FactorisedDistributionSampler sampler(distribution);
        return makeChainStatsPair(sampler, nSamples);
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
                    for(int type = 0; type < AGENT::typeDomainSize; ++type) {
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
    void load(Archive &ar, const unsigned int version) {
        std::chrono::steady_clock::duration::rep durationCount;
        ar >> meanVariance >> variogram >> varioStride >> meanEndState >> durationCount;
        sampleTime = std::chrono::steady_clock::duration(durationCount);
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar << meanVariance << variogram << varioStride << meanEndState << sampleTime.count();
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

};



#endif //GLPKTEST_CHAINSTATS_H
