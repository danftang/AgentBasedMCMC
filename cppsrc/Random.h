//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_RANDOM_H
#define GLPKTEST_RANDOM_H

#include <random>
#include <iostream>
#include <cassert>

class Random {
public:
    static std::mt19937 gen;

    static double nextDouble(double start = 0.0, double end=1.0) {
        return std::uniform_real_distribution<double>(start,end)(gen);
    }

    static int nextInt(int until) {
        return nextInt(0,until);
    }

    static int nextInt(int from, int until) {
        assert(from < until);
        return std::uniform_int_distribution<int>(from, until-1)(gen);
    }

    static int nextIntFromDiscrete(const std::vector<double> &probabilities) {
        return nextIntFromDiscrete(probabilities.begin(), probabilities.end());
    }

    template<typename InputIterator>
    static int nextIntFromDiscrete(InputIterator begin, InputIterator end) {
        return std::discrete_distribution<int>(begin, end)(gen);
    }

    static int nextPoisson(double lambda) {
        return std::poisson_distribution(lambda)(gen);
    }

    static int nextBinomial(int nTirals, double p) {
        return std::binomial_distribution<int>(nTirals, p)(gen);
    }

};


#endif //GLPKTEST_RANDOM_H
