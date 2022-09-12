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
#include "SparseFunction.h"

template<class DOMAIN>
class FactorisedDistribution {
public:
    typedef SparseFunction<std::pair<double,bool>,const DOMAIN &> function_type;
    typedef DOMAIN domain_type;

    std::vector<function_type>          logFactors;
    int                                 domainDimension;

    virtual std::function<const DOMAIN &()> sampler() {
        // TODO: implement this
        return nullptr;
    }

    void addFactor(SparseFunction<std::pair<double,bool>,const DOMAIN &> factor) {
        logFactors.push_back(std::move(factor));
    }

    FactorisedDistribution<DOMAIN> &operator *=(const FactorisedDistribution<DOMAIN> &other) {
        logFactors.reserve(logFactors.size() + other.logFactors.size());
        logFactors.insert(logFactors.end(), other.logFactors.begin(), other.logFactors.end());
        if(other.domainDimension > domainDimension) domainDimension = other.domainDimension;
        return *this;
    }

    FactorisedDistribution<DOMAIN> operator *(const FactorisedDistribution<DOMAIN> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    FactorisedDistribution<DOMAIN> operator *(const FactorisedDistribution<DOMAIN> &factoredDist) const & {
        FactorisedDistribution<DOMAIN> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }

    FactorisedDistribution<DOMAIN> operator *(FactorisedDistribution<DOMAIN> &&factoredDist) const & {
        factoredDist *= *this;
        return std::move(factoredDist);
    }

    // returns probability at point X
//    double P(const std::vector<T> &X) const {
//        return exp(logP(X));
//    }

    double exactFactorValue(int factorIndex, const DOMAIN &X) const {
        std::pair<double,bool> factorVal = logFactors[factorIndex](X);
        return factorVal.second?factorVal.first:-std::numeric_limits<double>::infinity();
    }

    double widenedFactorValue(int factorIndex, const DOMAIN &X) const {
        return logFactors[factorIndex](X).first;
    }


    double logPexact(const DOMAIN &X) const {
        double logP = 0.0;
        for(int i=0; i < logFactors.size(); ++i) logP += exactFactorValue(i,X);
        return logP;
    }

    double logPwidened(const DOMAIN &X) const {
        double logP = 0.0;
        for(int i=0; i < logFactors.size(); ++i) logP += widenedFactorValue(i,X);
        return logP;
    }

};


#endif //ABMCMC_FACTORISEDDISTRIBUTION_H
