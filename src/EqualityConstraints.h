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

    EqualityConstraints<COEFF> &operator +=(const EqualityConstraints<COEFF> &other) {
        this->insert(this->end(), other.begin(), other.end());
        return *this;
    }

    int maxCoefficientIndex() const {
        int max = 0;
        for(const EqualityConstraint<COEFF> &constraint: *this) {
            max = std::max(max, constraint.maxCoefficientIndex());
        }
        return max;
    }


    template<typename T>
    bool isValidSolution(const T &X) const {
        for(int i=0; i < this->size(); ++i) {
            if(!(*this)[i].isValidSolution(X)) return false;
        }
        return true;
    }
};

template<typename COEFF>
std::ostream &operator <<(std::ostream &out, const EqualityConstraints<COEFF> &constraints) {
    for(const EqualityConstraint<COEFF> &constraint: constraints) out << constraint << std::endl;
    return out;
}

#endif //ABMCMC_EQUALITYCONSTRAINTS_H
