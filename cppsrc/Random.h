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

};


#endif //GLPKTEST_RANDOM_H
