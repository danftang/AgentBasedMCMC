// A distribution defined inside a convex polyhedron
// L <= CX <= H
// with a probability that can be factorized into functions
// of the auxiliary variables. i.e. expressed in the form
//
// P(X) = \prod_i F_i(C_iX)
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_FACTOREDCONVEXDISTRIBUTION_H
#define ABMCMC_FACTOREDCONVEXDISTRIBUTION_H

#include <functional>
#include <cmath>
#include "version1/ConvexPolyhedron.h"

template<class T>
class FactoredConvexDistribution {
public:
    ConvexPolyhedron<T>                     constraints;
    std::vector<std::function<double(T)>>   factors;

    void addFactor(Constraint<T> constraint, std::function<double(T)> factor = nullLogProb) {
        constraints.push_back(std::move(constraint));
        factors.push_back(std::move(factor));
    }

    void reserve(size_t size) {
        constraints.reserve(size);
        factors.reserve(size);
    }

    FactoredConvexDistribution<T> &operator *=(const FactoredConvexDistribution<T> &other) {
        constraints += other.constraints;
        factors.reserve(factors.size() + other.factors.size());
        factors.insert(factors.end(), other.factors.begin(), other.factors.end());
        return *this;
    }

    FactoredConvexDistribution<T> operator *(const FactoredConvexDistribution<T> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    FactoredConvexDistribution<T> operator *(const FactoredConvexDistribution<T> &factoredDist) const & {
        FactoredConvexDistribution<T> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }

    FactoredConvexDistribution<T> operator *(FactoredConvexDistribution<T> &&factoredDist) const & {
        factoredDist *= *this;
        return std::move(factoredDist);
    }

    // returns probability at point X
    double P(const std::vector<T> &X) const {
        return exp(logP(X));
    }


    double logP(const std::vector<T> &X) const {
        double logP = 0.0;
        for(int i=0; i < factors.size(); ++i) {
            T Ai = constraints[i].coefficients * X;
            if(constraints[i].nObserved > Ai || Ai > constraints[i].upperBound) return -std::numeric_limits<double>::infinity();
            logP += factors[i](Ai);
        }
        return logP;
    }

protected:

    static double nullLogProb(T x) { return 0.0; }

};


template<class T>
std::ostream &operator <<(std::ostream &out, const FactoredConvexDistribution<T> &distribution) {
    for(int i=0; i < distribution.constraints.size(); ++i) {
        const Constraint<T> &constraint = distribution.constraints[i];
        out << constraint << "\t["
            << exp(distribution.factors[i](constraint.lowerBound)) << " ... "
            << exp(distribution.factors[i](constraint.upperBound)) << "]" << std::endl;
    }
    return out;
}

#endif //ABMCMC_FACTOREDCONVEXDISTRIBUTION_H
