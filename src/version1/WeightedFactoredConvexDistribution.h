// A distribution defined inside a convex polyhedron
// L <= CX <= H
// with a probability of the form
//
// P(X) = P_i(X) W(X)
//
// Where P_i(X) is factorized into functions of the auxiliary
// variables. i.e. expressed in the form
//
// P_i(X) = \prod_i F_i(C_iX)
//
// and W(X) is an importance weight function that can be updated
// by perturbations.
//
// IMPORTANCEFUNC should be a class which is a perturbable distribution with
// logImportance() and perturb(SparseVec) member functions. The template acts as
// a tag, the class will be instantiated once we convert to a sampler.
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_WEIGHTEDFACTOREDCONVEXDISTRIBUTION_H
#define ABMCMC_WEIGHTEDFACTOREDCONVEXDISTRIBUTION_H

#include "../FactoredConvexDistribution.h"
#include "PerturbableFunction.h"

template<class T>
class WeightedFactoredConvexDistribution: public FactoredConvexDistribution<T> {
public:

    std::function<std::unique_ptr<PerturbableFunction<T,double>>()> perturbableFunctionFactory;

    WeightedFactoredConvexDistribution()=default;

    explicit WeightedFactoredConvexDistribution(std::function<std::unique_ptr<PerturbableFunction<T,double>>()> perturbableFunctionFactory):
    perturbableFunctionFactory(perturbableFunctionFactory) {}

    WeightedFactoredConvexDistribution<T> operator *(const FactoredConvexDistribution<T> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    WeightedFactoredConvexDistribution<T> operator *(const FactoredConvexDistribution<T> &factoredDist) const & {
        WeightedFactoredConvexDistribution<T> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return std::move(copyOfThis);
    }


    WeightedFactoredConvexDistribution<T> &operator *=(const FactoredConvexDistribution<T> &factoredDist) {
        FactoredConvexDistribution<T>::operator*=(factoredDist);
        return *this;
    }

    friend WeightedFactoredConvexDistribution<T> operator *(const FactoredConvexDistribution<T> &lhs, WeightedFactoredConvexDistribution<T> &&rhs) {
        return rhs * lhs;
    }

    friend WeightedFactoredConvexDistribution<T> operator *(const FactoredConvexDistribution<T> &lhs, const WeightedFactoredConvexDistribution<T> &rhs) {
        return rhs * lhs;
    }

    double P(const std::vector<T> &X) const {
        return exp(logP(X));
    }

    double logP(const std::vector<T> &X) const {
        return FactoredConvexDistribution<T>::logP(X) + logImportance(X);
    }

    double logImportance(const std::vector<T> &X) const {
        std::unique_ptr<PerturbableFunction<T,double>> importance = perturbableFunctionFactory();
        importance->setState(X);
        return(importance->getLogValue(X));
    }

};


#endif //ABMCMC_WEIGHTEDFACTOREDCONVEXDISTRIBUTION_H
