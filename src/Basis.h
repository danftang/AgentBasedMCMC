// Represents a basis on a domain
// the basis may or may not span the domain
//
// Created by daniel on 26/09/22.
//

#ifndef ABMCMC_BASIS_H
#define ABMCMC_BASIS_H

#include <vector>
#include <assert.h>
#include <boost/serialization/access.hpp>

#include "include/SparseVec.h"
#include "include/debug.h"
#include "ConstrainedFactorisedDistribution.h"
#include "TableauEntropyMaximiser.h"
#include "TableauNormMinimiser.h"

template<class DOMAIN>
class Basis {
public:
    typedef typename DOMAIN::value_type value_type;

    std::vector<SparseVec<value_type>>  basisVectors;
    DOMAIN                              origin;

    Basis()=default;
//    Basis(const Basis<DOMAIN> &)=default;
//    Basis(Basis<DOMAIN> &&)=default;

    Basis(std::vector<SparseVec<value_type>> basisVecs, DOMAIN origin) :
            basisVectors(std::move(basisVecs)),
            origin(std::move(origin)) {
    }


    Basis(const ConstrainedFactorisedDistribution<DOMAIN> &distribution) {
        TableauNormMinimiser<value_type> tableau(distribution);
        basisVectors = tableau.getBasisVectors();
        origin += tableau.getOrigin(); // assumes default constructor of DOMAIN is the zero vector
        debug(sanityCheck(distribution.constraints));
    }

    Basis(const std::vector<double> &entropiesByVarIndex, const ConstrainedFactorisedDistribution<DOMAIN> &distribution, double constraintKappa) {
        TableauEntropyMaximiser<value_type> tableau(entropiesByVarIndex, distribution.constraints, constraintKappa);
        tableau.factorise();
        basisVectors = tableau.getBasisVectors();
        origin += tableau.getOrigin(); // assumes default constructor of DOMAIN is the zero vector
        debug(sanityCheck(distribution.constraints));
        debug(std::cout << "Created basis with entropy " << calculateEntropy(entropiesByVarIndex) << std::endl);
    }


    DOMAIN basisToDomain(const std::vector<value_type> &basisCoordinates) const {
        assert(basisCoordinates.size() == basisVectors.size());
        DOMAIN result(origin);
        for(int i=0; i<basisVectors.size(); ++i) result += basisVectors[i] * basisCoordinates[i];
        return result;
    }

    DOMAIN basisToDomain(std::function<value_type(int)> basisIndexToCoordinate) const {
        DOMAIN result(origin);
        for(int i=0; i<basisVectors.size(); ++i) result+= basisVectors[i] * basisIndexToCoordinate(i);
        return result;
    }

    // Calculates the expected entropy of an MCMC process on a fully factorised distribution
    // constrained to this basis, where the transition probability P(a,b) propto sqrt(P(b)/P(a))
    // entropiesByVar give the expected entropies of a unit perturbation in each dimension, defined as
    // E(i) = sum_x sqrt(P(x)P(x+1_i))
    // where 1_i is a vector with a single 1 in the i'th index, 0 elsewhere
    double calculateEntropy(const std::vector<double> &entropiesByVar) {
        double S = 0.0;
        for(const auto &basisVector: basisVectors) {
            double basisEntropy = 1.0;
            for(const auto &sparseEntry :basisVector) {
                basisEntropy *= pow(entropiesByVar[sparseEntry.index()], abs(sparseEntry.value()));
            }
            S += basisEntropy;
        }
        return S;
    }

private:
    friend class boost::serialization::access;

    void sanityCheck(const EqualityConstraints<value_type> &constraints) {
        // test that the basis is truly a basis of the constraints in the distribution
        assert(origin.size() == DOMAIN::dimension);
        assert(basisVectors.size() == DOMAIN::dimension - constraints.size());
        assert(constraints.isValidSolution(origin));
        for (int i = 0; i < basisVectors.size(); ++i) {
            origin += basisVectors[i];
            assert(constraints.isValidSolution(origin));
            origin -= basisVectors[i];
        }
    }


    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & basisVectors & origin;
    }
};


#endif //ABMCMC_BASIS_H
