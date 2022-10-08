// Represents a linear equation
// MX = E
// where M is a matrix.
//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_EQUALITYCONSTRAINTS_H
#define ABMCMC_EQUALITYCONSTRAINTS_H

#include <vector>
#include "EqualityConstraint.h"

template<class COEFF>
class EqualityConstraints: public std::vector<EqualityConstraint<COEFF>> {
public:
    typedef COEFF coefficient_type;

    EqualityConstraints<COEFF> &operator+=(const EqualityConstraints<COEFF> &other) {
        this->insert(this->end(), other.begin(), other.end());
        return *this;
    }

    int maxCoefficientIndex() const {
        int max = 0;
        for (const EqualityConstraint<COEFF> &constraint: *this) {
            max = std::max(max, constraint.maxCoefficientIndex());
        }
        return max;
    }


    template<typename T>
    bool isValidSolution(const T &X) const {
        for (int i = 0; i < this->size(); ++i) {
            if (!(*this)[i].isValidSolution(X)) return false;
        }
        return true;
    }


    // transform the constraint coefficients to a new domain, given that there is a linear transform
    // from the new domain to the old domain.
    // A = MB, where linearTransform is a vector of rows of M
    // So that NA = (NM)B
    EqualityConstraints<COEFF> toDomain(const std::vector<SparseVec<COEFF>> &linearTransform) {
        EqualityConstraints<COEFF> newConstraints;
        for(auto &thisConstraint: *this) {
            std::map<int,COEFF> newCoefficients;
            for(const auto &entry: thisConstraint.coefficients) {
                for(const auto &transformEntry: linearTransform[entry.index()]) {
                    COEFF increment = entry.value() * transformEntry.value();
                    auto ins = newCoefficients.try_emplace(transformEntry.index(), increment);
                    if(!ins.second) ins.first += increment;
                }
            }
            newConstraints.emplace_back(SparseVec(newCoefficients), thisConstraint.constant);
        }
        return newConstraints;
    };
};

template<typename COEFF>
std::ostream &operator <<(std::ostream &out, const EqualityConstraints<COEFF> &constraints) {
    for(int i=0; i<constraints.size() && i<100; ++i) {
        out << constraints[i] << std::endl;
    }
    if(constraints.size() >100) out << "...  ..." << std::endl;
    return out;
}

#endif //ABMCMC_EQUALITYCONSTRAINTS_H
