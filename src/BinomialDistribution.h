//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_BINOMIALDISTRIBUTION_H
#define GLPKTEST_BINOMIALDISTRIBUTION_H

#include <boost/math/distributions/binomial.hpp>
#include "Random.h"
#include "Distribution.h"
#include "ConvexPolyhedron.h"
#include "ConvexPMF.h"
#include "SampleStatistics.h"
#include "debug.h"

class BinomialDistribution: public Distribution<std::vector<double>> {
public:
    std::vector<boost::math::binomial_distribution<double>> binomials;

    explicit BinomialDistribution(std::vector<boost::math::binomial_distribution<double>> binoms): binomials(std::move(binoms)) {}

    // set up binomials that have given mean and approximate variance.
    // mean, m = np
    // variance, v = np(1-p)
    // since n must be an integer, calculate n first
    // n = m/p = m/(1-v/m) = m^2/(m-v)
    // p = m/n
    // this gives Delta_v approx= Delta_n d(m-m^2/n)/dn = Delta_n -m^2 ln(n)
    // TODO: Use maximisation of likelihood/posterior?
    explicit BinomialDistribution(const SampleStatistics &stats): binomials(stats.nDimensions()) {
        std::valarray<double> means = stats.means();
//        std::valarray<double> variances = stats.variances();
        for(int i=0; i<stats.nDimensions(); ++i) {
//            double nExact = means[i]*means[i]/(means[i] - variances[i]);
            int n = std::max(1.0,std::round(stats.max[i] + 1.0));//std::max(1.0,std::ceil(nExact)); // std::max(1.0, std::round(nExact));
            double p = std::min(means[i]/n,1.0);
//            debug(std::cout << "Creating binomial (" << n << "," << p << ")" << std::endl);
            binomials[i] = boost::math::binomial_distribution<double>(n, p);
        }
    }

//    double p(int i) const { return weights.at(i); }
//    double &p(int i) { return weights[i]; }

    double logP(const std::vector<double> &X) const {
        double logP = 0.0;
        for(int i=0; i < nDimensions(); ++i) {
            double successes = fabs(X[i]);
            if(successes <= binomials[i].trials()) {
                //            logP +=  successes*log(p) + (trials-successes)*log(1.0-p) + lgamma(trials+1) - lgamma(trials-successes+1) - lgamma(successes+1); // fails when p=0.0 or p=1.0
                logP += log(boost::math::pdf(binomials[i], successes));
            } else {
                logP -= INFINITY;
                i = nDimensions();
            }
        }
//        std::cout << "Binomial logP of " << X << " = " << logP << "  P = " << exp(logP) << std::endl;
        return logP;
    }

    std::vector<double> nextSample() const {
        std::vector<double> state(binomials.size());
        for(int i=0; i<binomials.size(); ++i) {
            state[i] = Random::nextBinomial(binomials[i].trials(), binomials[i].success_fraction());
//            std::cout << "Occupation " << i <<  " is " << state[i] << std::endl;
        }
        return state;
    }

    ConvexPMF<std::vector<double>> PMF() const {
        return ConvexPMF<std::vector<double>>([*this](const std::vector<double> &X) {
            return logP(X);
        }, nDimensions(), convexSupport());
    }

    ConvexPolyhedron convexSupport() const {
        ConvexPolyhedron support;
        for(int d=0; d < nDimensions(); ++d) {
            Constraint constraint = (0.0 <= 1.0*X(d) <= binomials[d].trials());
            if(binomials[d].success_fraction() == 0.0) {
                constraint.upperBound = 0.0;
            } else if(binomials[d].success_fraction() == 1.0) {
                constraint.lowerBound = binomials[d].trials();
            }
            support.push_back(constraint);
        }
        return support;
    }

    std::function<std::vector<double>()> sampler() const {
        return [*this]() { return this->nextSample(); };
    }

    int nDimensions() const { return binomials.size(); }

    std::vector<double> means() const {
        std::vector<double> mu(nDimensions());
        for(int i=0; i < nDimensions(); ++i) {
            mu[i] = boost::math::mean(binomials[i]);
            assert(mu[i] == binomials[i].success_fraction() * binomials[i].trials());
        }
        return mu;
    }

    friend std::ostream &operator <<(std::ostream &out, const BinomialDistribution &binomial) {
        for(auto &binom : binomial.binomials) {
            out << "(" << binom.trials() << ", " << binom.success_fraction() << ") ";
        }
        out << std::endl;
        return out;
    }

};


#endif //GLPKTEST_BINOMIALDISTRIBUTION_H
