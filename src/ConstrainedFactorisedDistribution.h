// A ConstrainedFactorisedDistribution represents a probability distribution over
// the template class T of the form
//
// P(X) = \prod_i P_i(X)
//
// subject to
//
// MX = C
// X \in {set of integers}
//
// for some set of sparse, widened functions P_i.
// The template class T should have a subscript operator.
//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H
#define ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H

#include <vector>
#include <limits>
#include "EqualityConstraints.h"
#include "SparseWidenedFunction.h"
#include "FactorisedDistribution.h"

template<class T>
class ConstrainedFactorisedDistribution: public FactorisedDistribution<T> {
public:
    EqualityConstraints<T>                  constraints;        // linear constraints

    virtual std::function<const T &()> sampler() {
        // TODO: Implement this
        return nullptr;
    }

    void addConstraint(EqualityConstraint<T> constraint) {
        constraints.push_back(std::move(constraint));
    }

    ConstrainedFactorisedDistribution<T> &operator *=(const ConstrainedFactorisedDistribution<T> &other) {
        constraints += other.constraints;
        FactorisedDistribution<T>::operator *=(other);
        return *this;
    }

    ConstrainedFactorisedDistribution<T> &operator *=(const FactorisedDistribution<T> &other) {
        FactorisedDistribution<T>::operator *=(other);
        return *this;
    }


    ConstrainedFactorisedDistribution<T> operator *(const ConstrainedFactorisedDistribution<T> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    ConstrainedFactorisedDistribution<T> operator *(const FactorisedDistribution<T> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }


    ConstrainedFactorisedDistribution<T> operator *(const ConstrainedFactorisedDistribution<T> &factoredDist) const & {
        ConstrainedFactorisedDistribution<T> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }

    ConstrainedFactorisedDistribution<T> operator *(const FactorisedDistribution<T> &factoredDist) const & {
        ConstrainedFactorisedDistribution<T> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }


    ConstrainedFactorisedDistribution<T> operator *(ConstrainedFactorisedDistribution<T> &&factoredDist) const & {
        factoredDist *= *this;
        return std::move(factoredDist);
    }

    friend ConstrainedFactorisedDistribution<T> operator *(const FactorisedDistribution<T> &fDist, const ConstrainedFactorisedDistribution<T> &cfDist) {
        return cfDist * fDist;
    }

    friend ConstrainedFactorisedDistribution<T> operator *(const FactorisedDistribution<T> &fDist, ConstrainedFactorisedDistribution<T> &&cfDist) {
        return std::move(cfDist) * fDist;
    }


    // returns probability at point X
//    double P(const std::vector<T> &X) const {
//        return exp(logP(X));
//    }


    double logPexact(const std::vector<T> &X) const {
        double logP = 0.0;
        for(int i=0; i < this->logFactors.size(); ++i) {
            if(constraints[i].coefficients * X != constraints[i].constant) return -std::numeric_limits<double>::infinity();
            logP += this->logFactors[i].exactValue(X);
        }
        return logP;
    }

};


#endif //ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H
