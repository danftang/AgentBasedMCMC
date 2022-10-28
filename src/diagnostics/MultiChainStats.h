// MCMC convergence statistics
//
// Created by daniel on 04/11/2021.
//

#ifndef GLPKTEST_MULTICHAINSTATS_H
#define GLPKTEST_MULTICHAINSTATS_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <future>

#include "ChainStats.h"
#include "../include/StlStream.h"

class MultiChainStats: public std::vector<ChainStats> {
public:
//    std::string cpuinfo;
//    double kap0pa;
//    std::string problemDescription;

    MultiChainStats()=default;

    // Generates two chain stats from first and last half of a 2xhalfSize chain
    template<class SAMPLER>
    MultiChainStats(int nSamples, SAMPLER &sampler, int nLags = ChainStats::defaultNLags, double maxLagProportion = ChainStats::defaultMaxLagProportion) {
        push_back(ChainStats(nSamples/2, sampler, nLags, maxLagProportion));
        push_back(ChainStats(nSamples/2, sampler, nLags, maxLagProportion));
    }

    MultiChainStats(std::vector<MultiChainStats> otherStats) {
        for(MultiChainStats &mcs: otherStats) {
            for(ChainStats &chain: mcs) {
                push_back(std::move(chain));
            }
        }
    }
//    MultiChainStats(double Kappa, std::string description=""):
//        kappa(Kappa),
//        problemDescription(std::move(description)) {
//
//    }

//    MultiChainStats &operator +=(std::vector<ChainStats> &&chains) {
//        reserve(size()+chains.size());
//        for(ChainStats &chain: chains) {
//            if(this->size()>0) assert(chain.nSamples() == this->nSamplesPerChain()); // all chains should be of the same length
//            this->push_back(std::move(chain));
//        }
//        return *this;
//    }
//
//    MultiChainStats &operator +=(const std::vector<ChainStats> &chains) {
//        reserve(size()+chains.size());
//        for(const ChainStats &chain: chains) {
//            if(this->size()>0) assert(chain.nSamples() == this->nSamplesPerChain()); // all chains should be of the same length
//            this->push_back(chain);
//        }
//        return *this;
//    }
    MultiChainStats &add(std::pair<ChainStats,ChainStats> &&chains) {
        reserve(size()+2);
        if(size()>0) {
            assert(chains.first.nSamples() == nSamplesPerChain()); // all chains should be of the same length
            assert(chains.second.nSamples() == nSamplesPerChain());
        }
        push_back(std::move(chains.first));
        push_back(std::move(chains.second));
        return *this;
    }


    int nSamplesPerChain() const {
        return size()==0?0:front().nSamples();
    }

    int nChains() const { return size(); }

    int dataDimension() const {
        return size()==0?0: front().dataDimension();
    }

    std::valarray<double> meanEndState() const {
        std::valarray<double> mean = (*this)[0].meanEndState;
        for(int j=1; j<size(); ++j) {
            mean += (*this)[j].meanEndState;
        }
        mean /= size();
        return mean;
    }

    // Calculates the effective number of samples as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, page 287.
    //
    // n_eff = mn/(1 + 2Dt sum_t rho_t)
    // where
    // Dt is the stride of t in the sum
    // rho_t = 1 - V_t/2var+
    // where Vt is the mean variogram over all chains at time-lag t
    // and var+ is an overestimate of the variance of the sampled distribution as defined below
    // and the sum over t starts at 0 and continues until rho_t + rho_t-1 < 0 for some odd t.
    std::valarray<double> effectiveSamples() const {
        std::valarray<double> neff(dataDimension());
        auto V = meanVariogram();
        std::valarray<double> varplus = varPlus();
        // calculate one dimension at a time, as each dimension may have a different stopping t
        for(int d=0; d < dataDimension(); ++d) {
            double vp = varplus[d];
            int t = 0;
            double nextS = 1.0 - V[0][d]/(2.0*vp);
            double sumt = 0.0;
            while (nextS > 0.0) {
                sumt += nextS;
                t += 1;
                nextS = std::min(
                        nextS,
                        (t < V.size())?(1.0 - V[t][d]/(2.0*vp)):-1.0
                        );
            }
            neff[d] = nSamplesPerChain() * nChains() / (1.0 + 2.0 * front().varioStride * sumt);
        }
        return neff;
    }


    // Approximates the autocorrelation averaged over all chains as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, page 286.
    //
    // rho_t = 1 - V_t/2var+
    //
    // where Vt is the mean variogram over all chains at time-lag t
    // and var+ is an overestimate of the variance of the sampled distribution as defined below
    std::valarray<std::valarray<double>> autocorrelation() const {
        auto V = meanVariogram();
        std::valarray<double> oneover2varplus = 1.0 / (varPlus()*2.0);
        std::valarray<std::valarray<double>> rhoHat( V.size());
        for(int t=0; t<V.size(); ++t) {
            rhoHat[t] = 1.0 - V[t] * oneover2varplus;
        }
        return rhoHat;
    }


    // The mean of the variograms in each chain
    std::valarray<std::valarray<double>> meanVariogram() const {
        std::valarray<std::valarray<double>> V(std::valarray<double>(0.0, dataDimension()), front().variogram.size());
        for(const ChainStats &chain: *this) {
            V += chain.variogram;
        }
        for(std::valarray<double> &Vt: V) {
            Vt /= nChains();
        }
        return V;
    }

    // The potential scale reduction for each dimension of the sample, as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, page 285.
    //
    // R_hat = sqrt((n-1)/n + B/(nW))
    std::valarray<double> potentialScaleReduction() const {
        return sqrt((B() / (W()*(1.0 * nSamplesPerChain()))) + ((nSamplesPerChain() - 1.0) / nSamplesPerChain()));
    }


    // Returns an over-estimate of the variance of the sampled distribution, as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, equation (11.3), page 284.
    //
    // var+ = (n-1/n) W + (1/n) B
    std::valarray<double> varPlus() const {
        return W() * ((nSamplesPerChain() - 1.0) / nSamplesPerChain()) + B() * (1.0 / nSamplesPerChain());
    }


    // Between sequence variance (i.e. the variance of the chain-means), defined as
    // B = n/(m-1) sum_m (mu_m - mu)^2
    // where n is number of samples, m is number of chains,
    // mu_m is the mean of the m'th chain and mu is the mean over all chains
    std::valarray<double> B() const {
        std::valarray<double> variance(0.0, dataDimension());
        std::valarray<double> mu = mean();
        for(const ChainStats &chain: *this) {
            std::valarray<double> Dmu_m = chain.meanVariance.mean() - mu;
            variance += Dmu_m * Dmu_m;
        }
        variance *= nSamplesPerChain() / (nChains() - 1);
        return variance;
    }

    // Within sequence variance (i.e. the mean of the chain-variances), defined as
    // W = 1/m sum_m s2_m
    // where m is the number of chains and s2_m is the variance of the m'th chain
    std::valarray<double> W() const {
        std::valarray<double> variance(0.0, dataDimension());
        for(const ChainStats &chain: *this) {
            variance += chain.meanVariance.sampleVariance();
        }
        variance /= nChains();
        return variance;
    }

    // mean over all samples and all chains
    std::valarray<double> mean() const {
        std::valarray<double> mu(0.0, dataDimension());
        for(const ChainStats &chain: *this) {
            mu += chain.meanVariance.mean();
        }
        mu /= nChains();
        return mu;
    }

    std::chrono::steady_clock::duration totalSampleTime() const {
        std::chrono::steady_clock::duration t(0);
        for(const ChainStats &chain: *this) {
            t += chain.sampleTime;
        }
        return t;
    }

//    static std::string getCpuInfo() {
//        FILE *pipe = popen("lscpu","r");
//        std::stringstream strstr;
//        while(!feof(pipe)) strstr << static_cast<char>(fgetc(pipe));
//        fclose(pipe);
//        return strstr.str();
//    }

    friend std::ostream &operator <<(std::ostream &out, const MultiChainStats &multiChainStats) {
//        out << "MultiChainStats for " << multiChainStats.problemDescription << " " << multiChainStats.nChains() << " chains with " << multiChainStats.nSamplesPerChain() << " samples" << std::endl;
        out << "MultiChainStats for "  << multiChainStats.nChains() << " chains with " << multiChainStats.nSamplesPerChain() << " samples" << std::endl;
//        out << multiChainStats.cpuinfo << std::endl;
        out << "Exec time = " << multiChainStats.totalSampleTime() << std::endl;
//        out << "W = " << multiChainStats.W() << std::endl;
//        out << "B = " << multiChainStats.B() << std::endl;
        out << "Samples   Sums    SumOfSquares    VarioStride     Vario" << std::endl;
        for(const ChainStats &chain: multiChainStats) {
            out << " " << chain.meanVariance.nSamples << " " << chain.meanVariance.sum << " " << chain.meanVariance.sumOfSquares << " " << chain.varioStride << " " << chain.variogram << std::endl;
//            std::cout << "MCMC stats:" << std::endl;
//            std::cout << chain.samplerStats << std::endl;
        }
        return out;
    }


//    template<class DISTRIBUTION>
//    static MultiChainStats analyseConvergence(const DISTRIBUTION &distribution, int nSamples, int nThreads = 4) {
//        MultiChainStats allThreadResults;
//        allThreadResults.reserve(2 * nThreads);
//
//        std::future<std::pair<ChainStats,ChainStats>> futureResults[nThreads-1];
//        for(int thread = 0; thread < nThreads-1; ++thread) {
//            futureResults[thread] = std::async(std::launch::async, [&distribution, nSamples]() {
//                return ChainStats::makeChainStatsPair(distribution, nSamples);
//            });
//        }
//        allThreadResults.add(ChainStats::makeChainStatsPair(distribution, nSamples));
//
//        for(int thread=0; thread<nThreads-1; ++thread) {
//            futureResults[thread].wait();
//            allThreadResults.add(futureResults[thread].get());
//        }
//
//        return allThreadResults;
//    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & static_cast<std::vector<ChainStats> &>(*this);
    }
};

#endif //GLPKTEST_MULTICHAINSTATS_H
