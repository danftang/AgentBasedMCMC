//
// Created by daniel on 04/11/2021.
//

#ifndef GLPKTEST_MULTICHAINSTATS_H
#define GLPKTEST_MULTICHAINSTATS_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "ChainStats.h"

class MultiChainStats: public std::vector<ChainStats> {
public:

    MultiChainStats &operator +=(MultiChainStats &&other) {
        reserve(size()+other.size());
        for(ChainStats &chain: other) {
            if(this->size()>0) assert(chain.nSamples() == this->nSamples()); // all chains should be of the same length
            this->push_back(std::move(chain));
        }
        return *this;
    }

    int nSamples() const {
        return size()==0?0:front().nSamples();
    }

    int nChains() const { return size(); }

    int dimension() const {
        return size()==0?0:front().dimension();
    }

    std::vector<double> meanEndState() {
        std::vector<double> mean = (*this)[0].meanEndState;
        for(int j=1; j<size(); ++j) {
            for(int i=0; i<mean.size(); ++i) {
                mean[i] += (*this)[j].meanEndState[i];
            }
        }
        for(int i=0; i<mean.size(); ++i) mean[i] /= size();
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
    std::valarray<double> effectiveSamples() {
        std::valarray<double> neff(dimension());
        auto V = meanVariogram();
        std::valarray<double> varplus = varPlus();
        // calculate one dimension at a time, as each dimension may have a different stopping t
        for(int d=0; d<dimension(); ++d) {
            double vp = varplus[d];
//            int t = 1;
//            double nextS = 2.0 - (V[1][d] + V[0][d])/(2.0*vp);
//            double sumt = 0.0;
//            while (nextS > 0.0) {
//                sumt += nextS;
//                t += 2;
//                nextS =(t < V.size())?(2.0 - (V[t][d] + V[t-1][d])/(2.0*vp)):-1.0;
//            }
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
            neff[d] = nSamples()*nChains()/(1.0 + 2.0*front().varioStride*sumt);
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
    std::valarray<std::valarray<double>> autocorrelation() {
        auto V = meanVariogram();
        std::valarray<double> oneover2varplus = 1.0 / (varPlus()*2.0);
        std::valarray<std::valarray<double>> rhoHat( V.size());
        for(int t=0; t<V.size(); ++t) {
            rhoHat[t] = 1.0 - V[t] * oneover2varplus;
        }
        return rhoHat;
    }


    // The mean of the variograms in each chain
    std::valarray<std::valarray<double>> meanVariogram() {
        std::valarray<std::valarray<double>> V(std::valarray<double>(0.0,dimension()), front().vario.size());
        for(const ChainStats &chain: *this) {
            V += chain.vario;
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
    std::valarray<double> potentialScaleReduction() {
        return sqrt((B() / (W()*(1.0*nSamples()))) + ((nSamples() -1.0)/nSamples()));
    }


    // Returns an over-estimate of the variance of the sampled distribution, as defined in
    // Gelman, A., Carlin, J.B., Stern, H.S. and Rubin, D.B., 2021. Bayesian data analysis 3rd Edition.
    // Chapter 11, equation (11.3), page 284.
    //
    // var+ = (n-1/n) W + (1/n) B
    std::valarray<double> varPlus() {
        return W() * ((nSamples()-1.0)/nSamples()) + B() * (1.0/nSamples());
    }


    // Between sequence variance (i.e. the variance of the chain-means), defined as
    // B = n/(m-1) sum_m (mu_m - mu)^2
    // where n is number of samples, m is number of chains,
    // mu_m is the mean of the m'th chain and mu is the mean over all chains
    std::valarray<double> B() const {
        std::valarray<double> variance(0.0,dimension());
        std::valarray<double> mu = mean();
        for(const ChainStats &chain: *this) {
            std::valarray<double> Dmu_m = chain.meanVariance.mean() - mu;
            variance += Dmu_m * Dmu_m;
        }
        variance *= nSamples()/(nChains()-1);
        return variance;
    }

    // Within sequence variance (i.e. the mean of the chain-variances), defined as
    // W = 1/m sum_m s2_m
    // where m is the number of chains and s2_m is the variance of the m'th chain
    std::valarray<double> W() const {
        std::valarray<double> variance(0.0,dimension());
        for(const ChainStats &chain: *this) {
            variance += chain.meanVariance.sampleVariance();
        }
        variance /= nChains();
        return variance;
    }

    // mean over all samples and all chains
    std::valarray<double> mean() const {
        std::valarray<double> mu(0.0,dimension());
        for(const ChainStats &chain: *this) {
            mu += chain.meanVariance.mean();
        }
        mu /= nChains();
        return mu;
    }

    friend std::ostream &operator <<(std::ostream &out, const MultiChainStats &multiChainStats) {
        out << "MultiChainStats of " << multiChainStats.nChains() << " chains with " << multiChainStats.nSamples() << " samples" << std::endl;
        out << "W = " << multiChainStats.W() << std::endl;
        out << "B = " << multiChainStats.B() << std::endl;
        out << "Samples   Sums    SumOfSquares    VarioStride     Vario" << std::endl;
        for(const ChainStats &chain: multiChainStats) {
            out << chain.meanVariance.nSamples << " " << chain.meanVariance.sum << " " << chain.meanVariance.sumOfSquares << " " << chain.varioStride << " " << chain.vario << std::endl;
        }
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & static_cast<std::vector<ChainStats> &>(*this);
    }
};

#endif //GLPKTEST_MULTICHAINSTATS_H
