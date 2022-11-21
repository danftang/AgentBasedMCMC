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
        debug(isCanonical());
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


    // ensures that the first non-zero entry in each basis is non-basic (i.e. that no other basis-vector
    // has a non-zero in this index.
    // returns true if a canonical form was possible, falso otherwise
    bool makeCanonical() {
        std::vector<int> basisIndex(DOMAIN::size(),-1); // basisVector that contains a given non-basic
        std::vector<int> entryIndex(DOMAIN::size(),-1); // the non-zero index of the entry that contains the non-basic
        for(int i=0; i<basisVectors.size(); ++i) {
            for(int nnz = 0; nnz < basisVectors[i].sparseSize(); ++nnz) {
                int domainIndex = basisVectors[i].indices[nnz];
                if(basisIndex[domainIndex] == -1) {
                    basisIndex[domainIndex] = i;
                    entryIndex[domainIndex] = nnz;
                } else {
                    basisIndex[domainIndex] = -2;
                }
            }
        }
        int nNonBasic = 0;
        for(int j=0; j<DOMAIN::size(); ++j) {
            if(basisIndex[j] >= 0 && entryIndex[j] != 0) {
                const SparseVec<value_type> &basisVec = basisVectors[basisIndex[j]];
                std::swap(basisVec.indices[entryIndex[j]], basisVec.indices[0]); // bring non-basic to front
                std::swap(basisVec.values[entryIndex[j]], basisVec.values[0]);
                ++nNonBasic;
            }
        }
        return nNonBasic >= basisVectors.size();
    }


//    // returns the highest index in any basis vector
//    int maxNonZeroIndex() {
//        int max = 0;
//        for(const auto &vec: basisVectors) max = std::max(max, vec.maxNonZeroIndex());
//        return max;
//    }

    bool isCanonical() const {
        std::vector<int> nonBasicVars(basisVectors.size());
        std::transform(basisVectors.begin(), basisVectors.end(), nonBasicVars.begin(), [](const SparseVec<value_type> &basisVec) {
            return basisVec.indices[0];
        });
        std::sort(nonBasicVars.begin(), nonBasicVars.end());
        return std::adjacent_find(nonBasicVars.begin(), nonBasicVars.end()) == nonBasicVars.end();
    }

private:
    friend class boost::serialization::access;

    void sanityCheck(const EqualityConstraints<value_type> &constraints) {
        // test that the basis is truly a basis of the constraints in the distribution

//        std::cout << "Checking basis. Basis vectors = " << std::endl;
//        std::cout << basisVectors << std::endl;
//        std::cout << "origin = " << origin << std::endl;
//        std::cout << "against constraints: " << std::endl;
//        std::cout << constraints << std::endl;

        assert(origin.size() == DOMAIN::dimension);
        assert(basisVectors.size() == DOMAIN::dimension - constraints.size());
        assert(constraints.isValidSolution(origin));
        for (int i = 0; i < basisVectors.size(); ++i) {
            origin += basisVectors[i];
//            std::cout << "Cheking basis " << basisVectors[i] << std::endl;
//            std::cout << "X = " << origin << std::endl;
            assert(constraints.isValidSolution(origin));
            origin -= basisVectors[i];
        }
        assert(isCanonical());
    }


    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & basisVectors & origin;
    }
};


#endif //ABMCMC_BASIS_H
