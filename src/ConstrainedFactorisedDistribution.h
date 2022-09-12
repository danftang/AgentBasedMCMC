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
#include "SparseFunction.h"
#include "FactorisedDistribution.h"

template<typename DOMAIN, typename CONSTRAINTCOEFF = typename subscript_operator_traits<DOMAIN>::base_type>
class ConstrainedFactorisedDistribution: public FactorisedDistribution<DOMAIN> {
public:
    typedef CONSTRAINTCOEFF coefficient_type;

    EqualityConstraints<CONSTRAINTCOEFF>                  constraints;        // linear constraints

    virtual std::function<const DOMAIN &()> sampler() {
        // TODO: Implement this
        return nullptr;
    }

    void addConstraint(EqualityConstraint<CONSTRAINTCOEFF> constraint) {
        this->domainDimension = std::max(this->domainDimension, constraint.maxCoefficientIndex());
        constraints.push_back(std::move(constraint));
    }

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &operator *=(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &other) {
        constraints += other.constraints;
        FactorisedDistribution<DOMAIN>::operator *=(other);
        return *this;
    }

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &operator *=(const FactorisedDistribution<DOMAIN> &other) {
        FactorisedDistribution<DOMAIN>::operator *=(other);
        return *this;
    }


    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const FactorisedDistribution<DOMAIN> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }


    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &factoredDist) const & {
        ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const FactorisedDistribution<DOMAIN> &factoredDist) const & {
        ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }


    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &&factoredDist) const & {
        factoredDist *= *this;
        return std::move(factoredDist);
    }

    friend ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const FactorisedDistribution<DOMAIN> &fDist, const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &cfDist) {
        return cfDist * fDist;
    }

    friend ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const FactorisedDistribution<DOMAIN> &fDist, ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &&cfDist) {
        return std::move(cfDist) * fDist;
    }


    // returns probability at point X
//    double P(const std::vector<T> &X) const {
//        return exp(logP(X));
//    }


    double logPexact(const DOMAIN &X) const {
        double logP = 0.0;
        if(!constraints.isValidSolution(X)) return -std::numeric_limits<double>::infinity();
        for(int i=0; i < this->logFactors.size(); ++i) {
            logP += this->exactFactorValue(i,X);
        }
        return logP;
    }

    double logPwidened(const DOMAIN &X) const {
        double logP = 0.0;
        double distanceToValidHyperplane = 0.0;
        for(const auto &constraint : constraints) {
            distanceToValidHyperplane += fabs(constraint.coefficients * X - constraint.constant);
        }
        logP += -ABM::kappa * distanceToValidHyperplane;
        for(int i=0; i < this->logFactors.size(); ++i) {
            logP += this->widenedFactorValue(i,X);
        }
        return logP;
    }

};

template<typename T>
std::ostream &operator <<(std::ostream &out, const ConstrainedFactorisedDistribution<T> &distribution) {
    out << distribution.constraints << std::endl;
    return out;
}

#endif //ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H
