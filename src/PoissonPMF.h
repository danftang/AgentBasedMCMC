//
// Created by daniel on 03/08/2021.
//

#ifndef GLPKTEST_POISSONPMF_H
#define GLPKTEST_POISSONPMF_H

#include <vector>
#include <cmath>
#include "Random.h"

// Represents a LogPMF where each variable is drawn from a Poisson distribution
class PoissonPMF {
public:
    std::vector<double> lambdas;
    std::vector<glp::Constraint> convexSupport;

    explicit PoissonPMF(std::vector<double> lambdas): lambdas(std::move(lambdas)) {}

    double lambda(int i) const { return lambdas.at(i); }
    double &lambda(int i) { return lambdas[i]; }

    double operator ()(const std::vector<double> &X) const {
        double logP = 0.0;
        for(int i=0; i < lambdas.size(); ++i) {
            double k = fabs(X[i]);
            double l = lambdas[i];
            logP += k*log(l) - l - lgamma(k+1); // log of Poisson
        }
        return logP;
    }

    std::vector<double> nextSample() {
        std::vector<double> state(lambdas.size());
        for(int i=0; i<lambdas.size(); ++i) {
            state[i] = Random::nextPoisson(lambda(i));
            std::cout << "Occupation " << i <<  " is " << state[i] << std::endl;
        }
        return state;
    }

};


#endif //GLPKTEST_POISSONPMF_H
