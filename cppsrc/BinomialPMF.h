//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_BINOMIALPMF_H
#define GLPKTEST_BINOMIALPMF_H

#include "Random.h"
#include "Distribution.h"

class BinomialPMF: public Distribution {
public:
    std::vector<double> weights;
    int trials;

    explicit BinomialPMF(std::vector<double> weights, int trials): weights(std::move(weights)), trials(trials) {}

    double p(int i) const { return weights.at(i); }
    double &p(int i) { return weights[i]; }

    double logP(const std::vector<double> &X) const {
        double logP = 0.0;
        for(int i=0; i < dimension(); ++i) {
            double successes = fabs(X[i]);
            if(successes <= trials) {
                double p = weights[i];
                //            logP +=  successes*log(p) + (trials-successes)*log(1.0-p) + lgamma(trials+1) - lgamma(trials-successes+1) - lgamma(successes+1); // fails when p=0.0 or p=1.0
                logP += log(boost::math::pdf(boost::math::binomial_distribution(1.0*trials, p), successes));
            } else {
                logP -= INFINITY;
                i = dimension();
            }
        }
//        std::cout << "Binomial logP of " << X << " = " << logP << "  P = " << exp(logP) << std::endl;
        return logP;
    }

    std::vector<double> nextSample() const {
        std::vector<double> state(weights.size());
        for(int i=0; i<weights.size(); ++i) {
            state[i] = Random::nextBinomial(trials,p(i));
//            std::cout << "Occupation " << i <<  " is " << state[i] << std::endl;
        }
        return state;
    }

    ConvexPMF PMF() const {
        return ConvexPMF([*this](const std::vector<double> &X) {
            return logP(X);
        },dimension(), convexSupport());
    }

    ConvexPolyhedron convexSupport() const {
        ConvexPolyhedron support;
        for(int d=0; d<dimension(); ++d) {
            glp::Constraint constraint = (0.0 <= 1.0*glp::X(d) <= trials);
            if(weights[d] == 0.0) {
                constraint.upperBound = 0.0;
            } else if(weights[d] == 1.0) {
                constraint.lowerBound = trials;
            }
            support.push_back(constraint);
        }
        return support;
    }

    std::function<std::vector<double>()> sampler() const {
        return [*this]() { return this->nextSample(); };
    }

    int dimension() const { return weights.size(); }

};


#endif //GLPKTEST_BINOMIALPMF_H
