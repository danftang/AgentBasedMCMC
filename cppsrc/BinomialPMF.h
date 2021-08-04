//
// Created by daniel on 04/08/2021.
//

#ifndef GLPKTEST_BINOMIALPMF_H
#define GLPKTEST_BINOMIALPMF_H

#include "Random.h"

class BinomialPMF {
public:
    std::vector<double> weights;
    std::vector<glp::Constraint> convexSupport;
    int trials;

    explicit BinomialPMF(std::vector<double> weights, int trials): weights(std::move(weights)), trials(trials) {}

    double p(int i) const { return weights.at(i); }
    double &p(int i) { return weights[i]; }

    double operator ()(const std::vector<double> &X) const {
        double logP = 0.0;
        for(int i=0; i < weights.size(); ++i) {
            double successes = fabs(X[i]);
            double p = weights[i];
            logP +=  successes*log(p) + (trials-successes)*log(1.0-p) + lgamma(trials+1) - lgamma(trials-successes+1) - lgamma(successes+1); // log of Binomial
        }
        return logP;
    }

    std::vector<double> nextSample() {
        std::vector<double> state(weights.size());
        for(int i=0; i<weights.size(); ++i) {
            state[i] = Random::nextBinomial(trials,p(i));
            std::cout << "Occupation " << i <<  " is " << state[i] << std::endl;
        }
        return state;
    }

};


#endif //GLPKTEST_BINOMIALPMF_H
