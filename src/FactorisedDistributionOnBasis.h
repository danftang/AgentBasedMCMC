// Represents a distribution over a linear subspace of DOMAIN. Each point on the
// subspace must be a linear combination of a set of basis vectors.
//
// The factors take DOMAIN co-ordinates, not basis coordinates
//
// Created by daniel on 02/10/22.
//

#ifndef ABMCMC_FACTORISEDDISTRIBUTIONONBASIS_H
#define ABMCMC_FACTORISEDDISTRIBUTIONONBASIS_H

#include "extratraits.h"
#include "SparseFunction.h"
#include "Basis.h"

template<typename DOMAIN, typename CONSTRAINTCOEFF = typename subscript_operator_traits<DOMAIN>::decay_type>
class FactorisedDistributionOnBasis {
public:
    typedef CONSTRAINTCOEFF coefficient_type;
    typedef SparseFunction<std::pair<double,bool>,const DOMAIN &> function_type;
    typedef DOMAIN domain_type;

    Basis<DOMAIN>               basis;
    std::vector<function_type>  factors;

    FactorisedDistributionOnBasis(Basis<DOMAIN> basis): basis(std::move(basis)) { }


    FactorisedDistributionOnBasis(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &distributionWithEqualityConstraints):
    basis(distributionWithEqualityConstraints),
    factors(distributionWithEqualityConstraints.factors) {}


    FactorisedDistributionOnBasis(ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &&distributionWithEqualityConstraints):
    basis(distributionWithEqualityConstraints),
    factors(std::move(distributionWithEqualityConstraints.factors)) { }

};


#endif //ABMCMC_FACTORISEDDISTRIBUTIONONBASIS_H
