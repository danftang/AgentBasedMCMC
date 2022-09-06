// Represents a probability distribution over template type T
// of the form:
//
// P(X) = \prod_i P_i(X)
//
// for some set of sparse, widened functions P_i.
// the template T should be a class with a subscript operator.
//
// Created by daniel on 02/09/22.
//

#ifndef ABMCMC_FACTORISEDDISTRIBUTION_H
#define ABMCMC_FACTORISEDDISTRIBUTION_H

#include <vector>
#include "SparseWidenedFunction.h"

template<class T>
class FactorisedDistribution {
public:
    std::vector<SparseWidenedFunction<double,T>>          logFactors;

    virtual std::function<const T &()> sampler() {
        // TODO: implement this
        return nullptr;
    }

    void addFactor(SparseWidenedFunction<double,T> factor) {
        logFactors.push_back(std::move(factor));
    }

    FactorisedDistribution<T> &operator *=(const FactorisedDistribution<T> &other) {
        logFactors.reserve(logFactors.size() + other.logFactors.size());
        logFactors.insert(logFactors.end(), other.logFactors.begin(), other.logFactors.end());
        return *this;
    }

    FactorisedDistribution<T> operator *(const FactorisedDistribution<T> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    FactorisedDistribution<T> operator *(const FactorisedDistribution<T> &factoredDist) const & {
        FactorisedDistribution<T> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }

    FactorisedDistribution<T> operator *(FactorisedDistribution<T> &&factoredDist) const & {
        factoredDist *= *this;
        return std::move(factoredDist);
    }

    // returns probability at point X
//    double P(const std::vector<T> &X) const {
//        return exp(logP(X));
//    }


    double logPexact(const T &X) const {
        double logP = 0.0;
        for(int i=0; i < logFactors.size(); ++i) logP += logFactors[i].exactValue(X);
        return logP;
    }

    double logPwidened(const T &X) const {
        double logP = 0.0;
        for(int i=0; i < logFactors.size(); ++i) logP += logFactors[i].widenedValue(X);
        return logP;
    }

};


#endif //ABMCMC_FACTORISEDDISTRIBUTION_H
