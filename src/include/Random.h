//
// Created by daniel on 11/05/2021.
//

#ifndef GLPKTEST_RANDOM_H
#define GLPKTEST_RANDOM_H

#include <random>
#include <iostream>
#include <cassert>
#include <chrono>
#include <mutex>

class Random {
protected:
    static std::mutex   seedMutex;
    static std::mt19937 seedGenerator;
public:
    thread_local static std::mt19937 gen;


    // thread-safe seed generator.
    // This can be used to uniquely seed different thread local instances of gen.
    static std::mt19937::result_type nextRandomSeed() {
        seedMutex.lock();
        std::mt19937::result_type seed = seedGenerator();
        seedMutex.unlock();
        return seed;
    }


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

    static bool nextBool() {
        return std::uniform_int_distribution<int>(0, 1)(gen) == 0;
    }

    static bool nextBool(double pTrue) {
        return std::bernoulli_distribution(pTrue)(gen);
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

inline std::mutex Random::seedMutex;
inline std::mt19937 Random::seedGenerator(
        static_cast<std::mt19937::result_type>(std::chrono::steady_clock::now().time_since_epoch().count()) +
        static_cast<std::mt19937::result_type>(reinterpret_cast<uintptr_t>(&seedMutex))); // attempt at random initialisation
inline thread_local std::mt19937 Random::gen;

#endif //GLPKTEST_RANDOM_H
